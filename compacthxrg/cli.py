"""The compactHxRG command"""
import argparse
import sys

from . import __version__
from .core import SLOPE_OFFSET, compact

EXAMPLES = """
examples:
  compactHxRG filename.fits                     writes filename_compact.fits
  compactHxRG filename.fits --overwrite True    replaces filename.fits
  compactHxRG *.fits.gz --outdir compact/       one compact file per input
"""


def str2bool(value):
    lower = value.lower()
    if lower in ('true', 't', 'yes', 'y', '1'):
        return True
    if lower in ('false', 'f', 'no', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError("expected True or False, got '%s'" % value)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='compactHxRG', epilog=EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Rewrite HxRG ramp-fit FITS files (extensions slope, intercept, '
                    'errors, n) about 5 times smaller, as small integers with '
                    'BSCALE/BZERO so that readers get the original values back.')
    parser.add_argument('files', nargs='+', metavar='filename.fits',
                        help='ramp-fit files, .fits or .fits.gz')
    parser.add_argument('--overwrite', type=str2bool, default=False, metavar='True/False',
                        help='True replaces each input file with its compact version; '
                             'False (default) writes <name>_compact.fits next to it')
    parser.add_argument('--suffix', default='_compact',
                        help="added to the file name when not overwriting "
                             "(default '_compact')")
    parser.add_argument('--outdir',
                        help='folder for the compact files when not overwriting '
                             '(default: next to each input)')
    parser.add_argument('--offset', type=int, default=SLOPE_OFFSET,
                        help='ADU added to slope * TINT so that values slightly below '
                             'zero are kept (default %d)' % SLOPE_OFFSET)
    parser.add_argument('--tint-key', default='TINT',
                        help='primary header keyword holding the integration time in '
                             's (default TINT)')
    parser.add_argument('--check', type=str2bool, default=True, metavar='True/False',
                        help='read each compact file back and print the quantization '
                             'error per extension (default True)')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args(argv)

    if args.overwrite and args.outdir:
        parser.error('--overwrite True replaces the input files, it cannot be '
                     'combined with --outdir')

    n_failed = 0
    for path in args.files:
        try:
            compact(path, overwrite=args.overwrite, suffix=args.suffix,
                    outdir=args.outdir, offset=args.offset, tint_key=args.tint_key,
                    run_check=args.check)
        except Exception as err:
            n_failed += 1
            print('%s: FAILED, %s' % (path, err), file=sys.stderr)
    if n_failed:
        print('%d of %d files failed' % (n_failed, len(args.files)), file=sys.stderr)
    return 1 if n_failed else 0


if __name__ == '__main__':
    sys.exit(main())
