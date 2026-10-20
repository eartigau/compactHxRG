# compactHxRG

`compactHxRG` rewrites HxRG ramp-fit FITS files about 5 times smaller by
storing every extension as small integers instead of 64-bit floats. `BSCALE`
and `BZERO` are set in each extension header, so any FITS reader (checked with
astropy and fitsio/CFITSIO) still gets the original values back, in the
original units, to within a rounding step that is well below the noise.

A 109 MB CPAPIR frame becomes 21 MB (8.3 MB if you also gzip it).

Une description du format en français est dans [FORMAT_fr.md](FORMAT_fr.md).

```
compactHxRG filename.fits
```

The compaction is lossy: values are rounded and out-of-range pixels are
clipped (details in [The format](#the-format)). Keep the originals somewhere
safe, and never run it with `--overwrite True` on archive data.

## Which files it works on

The input is a FITS file (`.fits` or `.fits.gz`) with:

- a primary HDU holding the header, with the integration time in seconds in a
  keyword (`TINT` by default, see `--tint-key`);
- four image extensions named `slope` (ADU/s), `intercept` (ADU), `errors`
  (ADU/s) and `n` (number of reads used in the fit).

This is the layout written by the HxRG ramp-fitting software, for example for
CPAPIR at the OMM 1.6 m from mid-2026 on. Files without these four extensions
are reported as failed and left alone.

## Step-by-step guide

### 1. Check your Python version

You need Python 3.9 or newer:

```
python3 --version
```

### 2. Install compactHxRG

Either install it directly from GitHub:

```
pip install git+https://github.com/eartigau/compactHxRG.git
```

or clone it first, if you want to read or modify the code (`-e` makes your
changes take effect without reinstalling):

```
git clone https://github.com/eartigau/compactHxRG.git
cd compactHxRG
pip install -e .
```

This installs `numpy` and `astropy` if they are missing, and creates the
`compactHxRG` command.

### 3. Check that the command works

```
compactHxRG --help
```

This should print the list of options. If you get `command not found`, the
folder where pip puts commands is not in your `PATH`; you can either add it,
or run the exact same thing as `python3 -m compacthxrg` (for example
`python3 -m compacthxrg --help`).

### 4. Compact one file and keep the original

```
compactHxRG filename.fits
```

This writes `filename_compact.fits` next to `filename.fits`, which is not
touched. The output looks like this:

```
filename.fits -> filename_compact.fits  (109.1 MB -> 21.0 MB, TINT = 2.621 s)
    slope     uint16 step     0.3815  max|err| 0.501 step  clipped 153  non-finite 3
    intercept uint8  step      50.95  max|err| 0.500 step  clipped 66565  non-finite 3
    errors    uint8  step     0.3815  max|err| 0.000 step  clipped 0  non-finite 0
    n         uint8  step          1  max|err| 0.000 step  clipped 0  non-finite 0
```

The compact file was read back and compared with the original, extension by
extension:

- `step`: the rounding step, in the units of the extension (ADU/s for `slope`
  and `errors`, ADU for `intercept`);
- `max|err|`: the largest difference with the original, in steps (0.5 is the
  expected maximum for rounding; the extra 0.001 is float32 precision when
  reading the file back);
- `clipped`: pixels outside the range that can be encoded, set to the nearest
  limit (hot pixels, dead pixels, pixels without valid reads);
- `non-finite`: NaN or inf pixels in the original, which read back as NaN.

### 5. Compact many files at once

Give as many files as you want, or a wildcard:

```
compactHxRG *.fits
compactHxRG *.fits.gz --outdir compact/
```

`--outdir` puts the compact files in another folder (created if needed)
instead of next to each input. If one file fails (not a ramp-fit file,
missing keyword...), the error is printed and the other files are still
processed.

### 6. Replace the original files (overwrite)

To save disk space on a local copy, replace each file by its compact version:

```
compactHxRG filename.fits --overwrite True
compactHxRG *.fits --overwrite True
```

The compact file is first written in a temporary folder and checked, and only
then takes the place of the original, so an interrupted run never leaves a
half-written file in place of the original. A `.fits.gz` input stays gzipped
(a 30.4 MB CPAPIR `.fits.gz` becomes 8.3 MB). Files that are already
compact are recognized and skipped, so running the command twice on the same
files does nothing the second time.

This cannot be undone: only use it on copies.

### 7. Read the compact files

Nothing special is needed; astropy applies `BSCALE`, `BZERO` and `BLANK` when
it reads the data:

```python
from astropy.io import fits

with fits.open('filename_compact.fits') as hdul:
    tint = hdul[0].header['TINT']
    slope = hdul['slope'].data          # ADU/s, NaN where the original was NaN/inf
    intercept = hdul['intercept'].data  # ADU
    errors = hdul['errors'].data        # ADU/s
    n = hdul['n'].data                  # number of reads

image_adu = slope * tint                # ADU
```

The values come back as float32 arrays instead of float64; float32 is far more
precise than the rounding step, so nothing is lost there.

Other readers give the same values, but not all of them turn `BLANK` pixels
into NaN: fitsio, for example, does not by default, and returns
`-1000 / TINT` ADU/s for the slope and one step above the top of the range for
the intercept at those pixels.

### 8. Update or uninstall

```
pip install --upgrade git+https://github.com/eartigau/compactHxRG.git
pip uninstall compactHxRG
```

(With a clone installed with `pip install -e .`, a `git pull` is enough to
update.)

## Options

| Option | Default | What it does |
|---|---|---|
| `--overwrite True/False` | `False` | `True` replaces each input file with its compact version; `False` writes `<name>_compact.fits` |
| `--suffix TEXT` | `_compact` | added to the file name when not overwriting |
| `--outdir FOLDER` | next to each input | where the compact files go when not overwriting (cannot be combined with `--overwrite True`) |
| `--offset ADU` | `1000` | added to `slope * TINT` before storing it, so that values down to `-offset + 1` ADU are kept (dark frames go slightly below zero) |
| `--tint-key KEYWORD` | `TINT` | primary header keyword holding the integration time in seconds |
| `--check True/False` | `True` | read each compact file back and print the comparison of step 4; `False` is a little faster |
| `--version` | | print the version |

The command returns exit code 0 if every file was compacted or skipped, 1 if
any failed, so it can be used in scripts.

From Python, the same thing is available as a function:

```python
from compacthxrg import compact

compact('filename.fits')                  # returns 'filename_compact.fits'
compact('filename.fits', overwrite=True)  # returns 'filename.fits'
```

## The format

The primary header is copied unchanged. Each extension keeps its name and
header; only the data type and the scaling keywords change.

| Extension | Original | Stored as | Rounding step | Range kept |
|---|---|---|---|---|
| `slope` | float64, ADU/s | uint16, `slope * TINT + 1000` | 1 ADU, that is `1 / TINT` ADU/s | -999 to 64535 ADU |
| `intercept` | float64, ADU | uint8, 0 to 254 | `2 (p99 - p1) / 254` ADU | `[p1, p99]` widened to twice its width |
| `errors` | float64, ADU/s | uint8, `errors * TINT` | 1 ADU, that is `1 / TINT` ADU/s | 0 to 255 ADU |
| `n` | uint16 | uint8 | exact | 0 to 255 |

Details:

- **slope**: `BSCALE = 1 / TINT` and `BZERO = (32768 - 1000) / TINT` (the
  32768 is how FITS stores unsigned 16-bit integers), so that the value read
  back is in ADU/s. Rounding to 1 ADU adds an error with an RMS of
  `1 / sqrt(12) = 0.29` ADU, well below the read noise. Code 0 is the `BLANK`
  value: NaN and inf pixels are stored as 0 and read back as NaN, so finite
  values below -999 ADU are clipped to -999 ADU.
- **intercept**: the 1st and 99th percentiles are measured on the pixels with
  `n > 0`, because pixels without any valid read hold a fill value
  (`3 * 65535` on CPAPIR) that can be more than 1% of the frame. The range from
  p1 to p99 is then doubled around its center and mapped onto codes 0 to 254.
  Code 255 is the `BLANK` value for NaN and inf. On CPAPIR frames the step is
  about 50 ADU, so the intercept is kept as a diagnostic rather than for
  precision work. (`BLANK = 0` would be valid FITS too, but astropy ignores a
  `BLANK` that is 0, hence 255.)
- **errors**: stored in ADU, like the slope, with `BSCALE = 1 / TINT`. Errors
  above 255 ADU, and NaN errors, are stored as 255.
- **n**: copied as integers, clipped at 255.

## Running the tests

```
pip install -e ".[test]"
pytest
```

The tests build small synthetic ramp-fit files, so no data is needed.
