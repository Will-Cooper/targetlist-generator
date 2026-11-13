# Targetlist Generator

Cut a target list down to what is observable on a given night, pair each target with a nearby telluric standard, and write the lists telescopes expect (IRTF starcat, ING object-visibility `.cat`) plus finder charts.
There is a point-and-click window for people who do not use Python, and a command line for those who do.

## Downloading (nothing else to install)
Download the file for your computer from the [latest release](https://github.com/Will-Cooper/targetlist-generator/releases/latest), unzip it, and open the app inside:
- **Windows**: `targetlist-gui-windows-x86_64.zip`, then `targetlist-gui\targetlist-gui.exe`.
- **macOS** (Apple Silicon): `targetlist-gui-macos-arm64.zip`, then `targetlist-gui.app`.
- **Linux**: `targetlist-gui-linux-x86_64.tar.gz`, then `targetlist-gui/targetlist-gui`.

The apps are not code-signed, so the first launch needs one extra step:
- **macOS** says the developer cannot be verified: open System Settings → Privacy & Security, scroll down and press **Open Anyway**, then confirm.
- **Windows** shows "Windows protected your PC": press **More info**, then **Run anyway**.

## Installing with uv (also gives the command line)
1. Install [uv](https://docs.astral.sh/uv/), a small tool that installs Python programs for you.
   - macOS / Linux, in a terminal: `curl -LsSf https://astral.sh/uv/install.sh | sh`
   - Windows, in PowerShell: `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`
2. Install the target list generator with its window:
   ```bash
   uv tool install "targetlist-generator[gui] @ git+https://github.com/Will-Cooper/targetlist-generator"
   ```
3. Start it with `targetlist-gui` (or `targetlist --help` for the command line).

Update later with `uv tool upgrade targetlist-generator`.
Python users can instead `pip install "targetlist-generator[gui] @ git+https://github.com/Will-Cooper/targetlist-generator"`.

## Using the window
Choose your target list, the evening the night starts and the observatory, then press **Run**.
The Targets tab shows what survived the cuts (with each target's best time, airmass, Moon distance and standard), the Airmass plot tab the night's airmass curves, and the Finder charts tab any charts made.
**Open output folder** takes you to the files.

## Using the command line
```bash
targetlist make targets.csv --date 2026-11-14 --site irtf
targetlist make targets.fits --date 2026-11-14 --lat 28.76 --lon -17.88 --max-airmass 1.8 --finders
targetlist make targets.csv -r 20h 4h -d -30 60 --no-standards
targetlist finder targets.csv J1023+2215 J0004-4044 --survey CDS/P/2MASS/J
```
`targetlist make -h` lists every option.

### Input
CSV, FITS, ECSV, VOTable, or a previous IRTF `.txt` output.
Column names are case-insensitive; required are `shortname` (or `name`), `ra`, `dec`, `spt` and a magnitude column (`--mag-col`, default `tmassk`).
RA and Dec can be decimal degrees or sexagesimal strings (`10:23:00`, `10h23m00s`); spectral types can be strings or numbers (M0 = 60, L0 = 70, T0 = 80, Y0 = 90).
`epoch`, `pmra` and `pmdec` are used when present (defaults 2000, 0, 0).

### Cuts
- **Night** (`--date` with `--site` or `--lat/--lon`): the night is when the Sun is below `--twilight` (default −12 deg).
  A target is kept if it spends at least `--min-time` minutes (default 30) at airmass ≤ `--max-airmass` (default 2).
  Site names come from astropy's registry, e.g. `irtf`, `lapalma`, `paranal`, `lasilla`, `keck`, `'Kitt Peak'`.
- **Moon** (with a night): `--min-moon-sep` drops targets closer than this many degrees to the Moon at mid-night (default off); the Moon's illumination is reported.
- **RA/Dec box** (`-r`, `-d`): degrees or e.g. `20h`; RA wraps through 0 when the lower limit is larger.

### Output (in `--out-dir`, default `targetlist_out/`)
- `<name>_irtf.txt`: IRTF starcat list, each target followed by its standard.
- `<name>_irtf_nostds.txt`: the targets only.
- `cats/<name>_cat_N.cat`: ING object-visibility lists, ten targets each.
- `fcharts/*.png`: finder charts (`--finders`, or `targetlist finder`).
- With a night: `<name>_schedule.csv` (each target's best UTC time, minimum airmass, first/last time below the airmass limit, Moon separation and standard) and `<name>_airmass.png` (airmass against UTC for every target).
  The lists are then ordered by best time; `--order input` keeps the input order.

### Standards
The bundled list holds 4969 B9V–A1V stars with K = 5–9 mag.
With a night, each target is paired with the standard nearest on the sky among those whose airmass `--std-offset` minutes (default 30) after the target's best time is within `--std-airmass-tol` (default 0.1) of the target's; if none is that close, the one closest in airmass is used.
Without a night, the standard nearest on the sky to where the target will be `--std-offset` minutes later in RA is used.
Choose the magnitude band with `--std-band`, a single spectral type with `--std-spt`, or your own list with `--standards-file`.

### Finder charts
Charts are cut from the CDS [hips2fits](https://alasky.cds.unistra.fr/hips-image-services/hips2fits) service (default survey 2MASS K, 2 arcmin), so no Java is needed.
Add another hips2fits service with `--hips2fits-url`.
If no hips2fits service answers and [Aladin Desktop](https://aladin.cds.unistra.fr/) is installed as an `aladin` command, the remaining charts are drawn in one Aladin session from the IRSA copies of the 2MASS, DSS2 and AllWISE surveys instead.

## Developing
```bash
pip install -e ".[dev,gui-test]"
pytest -m "not network"        # QT_QPA_PLATFORM=offscreen for the GUI tests on a headless machine
ruff check src tests
```
Version 1 (the single `targetlist_gen.py` script driving Aladin) is tagged `v1.0.4`.

## Citing
[![DOI](https://zenodo.org/badge/323075587.svg)](https://zenodo.org/badge/latestdoi/323075587)
