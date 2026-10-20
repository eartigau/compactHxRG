import os

import numpy as np
import pytest
from astropy.io import fits

from compacthxrg import compact
from compacthxrg.cli import main

TINT = 5.0


def make_ramp(path, shape=(64, 64), seed=0):
    """Small synthetic ramp-fit file with the CPAPIR layout"""
    rng = np.random.default_rng(seed)
    slope = rng.normal(100, 30, shape)  # ADU/s
    intercept = rng.normal(6000, 500, shape)
    errors = np.abs(rng.normal(2, 0.5, shape))
    n = np.full(shape, 15, dtype=np.uint16)
    n[3] = 0
    intercept[3] = 3 * 65535  # fill value of pixels without valid reads
    slope[0, 0], slope[0, 1] = np.nan, np.inf
    slope[1, 0], slope[1, 1] = 1e5, -1e3  # beyond the uint16 range once x TINT
    intercept[0, 0], intercept[0, 1] = np.nan, -np.inf
    errors[2, 2] = np.nan

    hdr = fits.Header()
    hdr['TINT'] = TINT
    hdr['OBJECT'] = 'test'
    fits.HDUList([fits.PrimaryHDU(header=hdr),
                  fits.ImageHDU(slope, name='slope'),
                  fits.ImageHDU(intercept, name='intercept'),
                  fits.ImageHDU(errors, name='errors'),
                  fits.ImageHDU(n, name='n')]).writeto(path)
    return dict(slope=slope, intercept=intercept, errors=errors, n=n)


def test_roundtrip(tmp_path):
    path = str(tmp_path / 'ramp.fits')
    orig = make_ramp(path)
    out = compact(path, verbose=False)
    assert out == str(tmp_path / 'ramp_compact.fits')
    assert os.path.exists(path)

    with fits.open(out) as hdul:
        assert hdul[0].header['OBJECT'] == 'test'
        bitpix = {name: hdul[name].header['BITPIX'] for name in orig}
        istep = hdul['intercept'].header['BSCALE']
        slope = hdul['slope'].data
        intercept = hdul['intercept'].data
        errors = hdul['errors'].data
        n = hdul['n'].data
    assert bitpix == dict(slope=16, intercept=8, errors=8, n=8)

    # half a step, plus the float32 rounding of the decoded values (BZERO is
    # about 31768 / TINT, so its float32 rounding is about 1e-3 of a step)
    tol = 0.505 / TINT
    good = np.isfinite(orig['slope'])
    good[1, :2] = False
    assert np.all(np.abs(slope[good] - orig['slope'][good]) <= tol)
    assert np.isnan(slope[0, 0]) and np.isnan(slope[0, 1])
    assert slope[1, 0] == pytest.approx(64535 / TINT)
    assert slope[1, 1] == pytest.approx(-999 / TINT)

    valid = (orig['n'] > 0) & np.isfinite(orig['intercept'])
    assert np.isnan(intercept[0, 0]) and np.isnan(intercept[0, 1])
    assert np.all(np.abs(intercept[valid] - orig['intercept'][valid]) <= 0.5 * istep + 1e-2)

    good = np.isfinite(orig['errors'])
    assert np.all(np.abs(errors[good] - orig['errors'][good]) <= tol)
    assert errors[2, 2] == pytest.approx(255 / TINT)
    assert np.array_equal(n, orig['n'])


def test_overwrite_then_skip(tmp_path):
    path = str(tmp_path / 'ramp.fits')
    make_ramp(path)
    assert compact(path, overwrite=True, verbose=False) == path
    with fits.open(path) as hdul:
        assert hdul['slope'].header['BITPIX'] == 16
    before = open(path, 'rb').read()
    assert compact(path, overwrite=True, verbose=False) is None
    assert open(path, 'rb').read() == before
    assert sorted(os.listdir(tmp_path)) == ['ramp.fits']


def test_gzip_input_and_outdir(tmp_path):
    path = str(tmp_path / 'ramp.fits.gz')
    make_ramp(path)
    out = compact(path, outdir=str(tmp_path / 'out'), verbose=False)
    assert out == str(tmp_path / 'out' / 'ramp_compact.fits')
    assert os.path.exists(out)


def test_refuses_to_replace_input_without_overwrite(tmp_path):
    path = str(tmp_path / 'ramp.fits')
    make_ramp(path)
    with pytest.raises(ValueError):
        compact(path, suffix='', verbose=False)


def test_cli(tmp_path, capsys):
    paths = [str(tmp_path / ('ramp%d.fits' % i)) for i in range(2)]
    for path in paths:
        make_ramp(path)
    assert main(paths + ['--overwrite', 'False']) == 0
    assert all(os.path.exists(p.replace('.fits', '_compact.fits')) for p in paths)
    assert main([str(tmp_path / 'missing.fits')]) == 1
    with pytest.raises(SystemExit):
        main(paths + ['--overwrite', 'True', '--outdir', str(tmp_path)])
    with pytest.raises(SystemExit):
        main(paths + ['--overwrite', 'maybe'])
