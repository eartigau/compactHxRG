"""
Compact HxRG ramp-fit files (extensions slope, intercept, errors, n) by storing
every extension as small integers, with BSCALE and BZERO set so that any FITS
reader gets back the original values (ADU/s for the slope and errors) to within
the quantization step.

    slope      uint16 = slope * TINT + 1000      1 ADU steps, -999 to 64535 ADU
    intercept  uint8, 0..254 spans [p1, p99] widened to twice its width
    errors     uint8 = errors * TINT             1 ADU steps, clipped at 255
    n          uint8

Non-finite slope/intercept pixels are stored as a BLANK code (0 for the slope,
255 for the intercept), read back as NaN; non-finite errors are stored as 255,
like errors too large to encode. The intercept percentiles leave out pixels
with n = 0, which hold a fill value. The primary header is copied as is.
"""
import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits

EXTENSIONS = ('slope', 'intercept', 'errors', 'n')
SLOPE_OFFSET = 1000  # ADU, room kept below zero for dark frames
UINT16_SHIFT = 32768  # FITS stores unsigned 16-bit as signed with BZERO = 32768
BITPIX_NAME = {8: 'uint8', 16: 'uint16'}  # as stored, before BSCALE/BZERO


def encode_slope(slope, tint, offset=SLOPE_OFFSET):
    """slope in ADU/s -> int16 raw values, BSCALE, BZERO, BLANK"""
    finite = np.isfinite(slope)
    code = np.zeros(slope.shape, dtype=np.int32)
    code[finite] = np.clip(np.rint(slope[finite] * tint + offset), 1, 65535)
    raw = (code - UINT16_SHIFT).astype(np.int16)
    # value = raw * BSCALE + BZERO = (code - offset) / tint
    return raw, 1.0 / tint, (UINT16_SHIFT - offset) / tint, -UINT16_SHIFT


def encode_intercept(intercept, n):
    """intercept -> uint8 codes 0..254 over [p1, p99] widened x2, 255 = BLANK

    (BLANK = 0 would be valid FITS, but astropy skips a BLANK that is 0)
    """
    finite = np.isfinite(intercept)
    # pixels without any valid read (n = 0) hold a fill value (3 * 65535) that
    # can be more than 1% of the frame, so they are left out of the percentiles
    valid = finite & (n > 0)
    if not np.any(valid):
        valid = finite
    if not np.any(valid):
        p1, p99 = 0.0, 0.0
    else:
        p1, p99 = np.percentile(intercept[valid], [1, 99])
    lo = p1 - (p99 - p1) / 2
    hi = p99 + (p99 - p1) / 2
    bscale = (hi - lo) / 254 if hi > lo else 1.0
    bzero = lo  # code 0 -> lo, code 254 -> hi
    code = np.full(intercept.shape, 255, dtype=np.uint8)
    code[finite] = np.clip(np.rint((intercept[finite] - bzero) / bscale), 0, 254)
    return code, bscale, bzero, 255


def encode_errors(errors, tint):
    """errors in ADU/s -> uint8 of errors * tint (ADU), clipped at 255"""
    finite = np.isfinite(errors)
    code = np.full(errors.shape, 255, dtype=np.uint8)
    code[finite] = np.clip(np.rint(errors[finite] * tint), 0, 255)
    return code, 1.0 / tint, 0.0, None


def make_hdu(orig_header, raw, bscale, bzero, blank):
    """Image HDU storing the raw integers with the given scaling keywords"""
    hdr = orig_header.copy()
    for key in ('BSCALE', 'BZERO', 'BLANK'):
        hdr.remove(key, ignore_missing=True)
    # astropy drops BSCALE/BZERO given with integer data; scale() is the way to
    # set them. It maps raw * bscale + bzero back onto raw exactly (rounding)
    hdu = fits.ImageHDU(data=raw * bscale + bzero, header=hdr)
    hdu.scale(raw.dtype.name, bscale=bscale, bzero=bzero)
    if blank is not None:
        hdu.header['BLANK'] = blank
    return hdu


def output_path(path, overwrite=False, suffix='_compact', outdir=None):
    """Where the compact version of path goes"""
    if overwrite:
        return path
    folder, base = os.path.split(path)
    for ext in ('.gz', '.fits'):
        if base.endswith(ext):
            base = base[:-len(ext)]
    return os.path.join(outdir or folder, base + suffix + '.fits')


def read_ramp(path, tint_key='TINT'):
    """Primary header, extension headers and data of a ramp-fit file

    Returns None for the data if the file is already compact.
    """
    with fits.open(path, memmap=False) as hdul:
        names = [hdu.name.lower() for hdu in hdul]
        missing = [name for name in EXTENSIONS if name not in names]
        if missing:
            raise ValueError('no %s extension, not an HxRG ramp-fit file'
                             % ', '.join(missing))
        if tint_key not in hdul[0].header:
            raise ValueError('no %s keyword in the primary header' % tint_key)
        primary = hdul[0].header.copy()
        headers = {name: hdul[name].header.copy() for name in EXTENSIONS}
        if headers['slope']['BITPIX'] > 0:
            return primary, headers, None
        data = {name: hdul[name].data for name in EXTENSIONS}
    return primary, headers, data


def check(orig, path, verbose=True):
    """Read the compact file back and compare it with the original values"""
    with fits.open(path) as hdul:
        for name in EXTENSIONS:
            # astropy removes BSCALE from the header once the data is read
            step = hdul[name].header.get('BSCALE', 1.0)
            dtype = BITPIX_NAME[hdul[name].header['BITPIX']]
            old = orig[name].astype(float)
            new = hdul[name].data.astype(float)
            finite = np.isfinite(old)
            if name in ('slope', 'intercept') and np.any(np.isnan(new) != ~finite):
                raise RuntimeError('%s: NaN pixels not restored' % name)
            diff = np.abs(new - old)
            with np.errstate(invalid='ignore'):  # NaN diffs, numpy < 1.19 warns
                ok = finite & (diff <= 0.505 * step)
            max_err = np.max(diff[ok]) / step if np.any(ok) else 0.0
            if verbose:
                print('    %-9s %-6s step %10.4g  max|err| %.3f step  clipped %d  '
                      'non-finite %d' % (name, dtype, step, max_err,
                                         np.sum(finite & ~ok), np.sum(~finite)))


def compact(path, overwrite=False, suffix='_compact', outdir=None,
            offset=SLOPE_OFFSET, tint_key='TINT', run_check=True, verbose=True):
    """Write the compact version of an HxRG ramp-fit file

    With overwrite=True the input file is replaced; otherwise the output is
    <name><suffix>.fits, next to the input or in outdir. Files that are already
    compact are skipped. Returns the output path, or None if skipped.
    """
    if not 0 <= offset <= 65535:
        raise ValueError('offset must be between 0 and 65535 ADU')
    out_path = output_path(path, overwrite, suffix, outdir)
    if not overwrite and os.path.abspath(out_path) == os.path.abspath(path):
        raise ValueError('the output name is the input name, use overwrite=True '
                         'to replace the input file')

    primary, headers, orig = read_ramp(path, tint_key)
    if orig is None:
        if verbose:
            print('%s is already compact, skipped' % path)
        return None
    tint = primary[tint_key]

    hdus = [fits.PrimaryHDU(header=primary)]
    hdus.append(make_hdu(headers['slope'], *encode_slope(orig['slope'], tint, offset)))
    hdus.append(make_hdu(headers['intercept'],
                         *encode_intercept(orig['intercept'], orig['n'])))
    hdus.append(make_hdu(headers['errors'], *encode_errors(orig['errors'], tint)))
    n = np.clip(orig['n'], 0, 255).astype(np.uint8)
    hdus.append(make_hdu(headers['n'], n, 1.0, 0.0, None))

    # written in a temporary folder under its final name (so a .gz stays
    # gzipped) and checked before it moves in place, possibly over the input
    folder = os.path.dirname(out_path) or '.'
    os.makedirs(folder, exist_ok=True)
    tmp_dir = tempfile.mkdtemp(prefix='.compactHxRG_', dir=folder)
    tmp_path = os.path.join(tmp_dir, os.path.basename(out_path))
    in_size = os.path.getsize(path)
    try:
        fits.HDUList(hdus).writeto(tmp_path)
        if verbose:
            print('%s -> %s  (%.1f MB -> %.1f MB, %s = %.3f s)'
                  % (path, out_path, in_size / 1e6, os.path.getsize(tmp_path) / 1e6,
                     tint_key, tint))
        if run_check:
            check(orig, tmp_path, verbose)
        os.replace(tmp_path, out_path)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    return out_path
