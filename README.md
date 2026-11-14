# Targetlist Generator

Plan an observing night from a target table.
Give it your targets (a file, or sent from TOPCAT), the date and the observatory, and it keeps what is visible, lays the targets through the night in blocks at low airmass, pairs each with a nearby telluric standard, and draws ING-style visibility plots.
It writes the target lists your telescope wants (IRTF, ING, Keck, Lick, ESO p2, Subaru, MMT, Magellan, or plain CSV/text/VOTable) and optional finder charts.
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
Work down the left-hand side:
1. **Target table**: browse to a file, or tick *Listen for TOPCAT (SAMP)* and use TOPCAT's *Interop → Send table to → Target List Generator* (a SAMP hub is started if none is running). The window checks the table has the columns it needs before you can plan.
2. **Time per target block**: the same observing time for every target (default 30 min) or a column of minutes from the table, plus the time for each standard (default 10 min) and for calibrations (default 5 min).
3. **Night and observatory**: pick the evening from the calendar and type part of an observatory's name or alias (e.g. *irtf*, *keck*, *palma*) or enter latitude/longitude; set the airmass limit and the Moon filter (automatic by default: 10 deg from a new Moon up to 40 deg from a full one).
4. **Telluric standards**: on or off, band, spectral type, airmass tolerance and maximum distance on the sky.
5. Optionally an **RA/Dec box**.
6. **Outputs**: the list format for the chosen observatory is ticked for you (otherwise CSV, or VOTable for a table sent from TOPCAT); tick others as you like, optionally finder charts, and choose where to save them (asked when you press *Plan the night* if not set).

A progress bar under *Plan the night* shows what is happening (finder charts can take a few minutes when drawn by Aladin).
The *Plan* tab shows the night: each target's block, its standard and how far the telescope moves to reach it.
The *Visibility* tab shows ING-style plots of every target's altitude through the night, with the planned blocks drawn thick, and the *Finder charts* tab each chart as it is made; both can be zoomed, panned and saved with the toolbar above them.

## Try it
`examples/example_targets.csv` holds 40 made-up targets (synthetic positions, not real objects) spread around the sky, with sexagesimal coordinates, string and numeric spectral types, proper motions and an `exptime` column of minutes.
Open it in the window (or send it from TOPCAT) and try, for example:
- the date picker and different observatories (*irtf*, *paranal*, *lapalma*) to see which targets survive and how the plan changes;
- *From column* → `exptime` instead of the same time for every target, or longer standard and calibration times, until the night overflows and some targets become `extra`;
- a stricter *Max airmass* or *Max distance* for standards, and the Moon filter set to *custom*;
- ticking extra output formats and finder charts.

From the command line: `targetlist make examples/example_targets.csv --date 2026-11-30 --site irtf --target-time-col exptime --formats irtf,keck,csv`.

## Using the command line
```bash
targetlist make targets.csv --date 2026-11-30 --site irtf
targetlist make targets.fits --date 2026-11-30 --site lapalma --target-time-col exptime --formats irtf,keck,votable
targetlist make targets.csv --date 2026-11-30 --lat 28.76 --lon -17.88 --min-moon-sep off --finders
targetlist make targets.csv -r 20h 4h -d -30 60 --no-standards
targetlist finder targets.csv J1023+2215 J0004-4044 --survey CDS/P/2MASS/J
targetlist sites keck
```
`targetlist make -h` lists every option.

### Input
CSV, FITS, ECSV, VOTable, or a previous IRTF `.txt` output.
Column names are case-insensitive; required are `shortname` (or `name`), `ra`, `dec`, `spt` and a magnitude column (`--mag-col`, default `tmassk`).
RA and Dec can be decimal degrees or sexagesimal strings (`10:23:00`, `10h23m00s`); spectral types can be strings or numbers (M0 = 60, L0 = 70, T0 = 80, Y0 = 90).
`epoch`, `pmra` and `pmdec` are used when present (defaults 2000, 0, 0); proper motions are mas/yr with `pmra` including cos(Dec), as in Gaia.

### Planning the night
With `--date` and an observatory (`--site`, or `--lat/--lon`), the night runs while the Sun is below `--twilight` (default −12 deg) and targets are kept if they spend at least `--min-time` minutes (default 30) at airmass ≤ `--max-airmass` (default 2) and are far enough from the Moon (`--min-moon-sep auto|off|<deg>`).
Each target then gets a block: its observing time (`--target-time`, or `--target-time-col`), calibrations (`--cal-time`) and its standard (`--std-time`).
Blocks are laid back to back through the night, each time taking the target closest to its best airmass of the night, so the list comes out in time (and so roughly RA) order.
No visible target is dropped: when the night is too full, the targets that do not fit are listed as `extra` at their best time, for you to choose between on the night.
Observatories come from a bundled copy of astropy's site registry (no internet needed); `targetlist sites <search>` lists them.

### Standards
The bundled list holds 4969 B9V–A1V stars with K = 5–9 mag.
Each standard is chosen for the time it will be observed, straight after its target: the nearest on the sky within `--std-max-sep` (default 20 deg) whose airmass is within `--std-airmass-tol` (default 0.1) of the target's.
Short slews come first: when no nearby standard matches the airmass, the closest match within `--std-max-sep` is used, and only when there is no standard that near at all is one further away taken; both compromises are noted in the plan's `std_note`.
Without a night, the standard nearest on the sky to where the target will be after its block is used.
Choose the magnitude band with `--std-band`, a single spectral type with `--std-spt`, or your own list with `--standards-file`.

### Output (in `--out-dir`, default `targetlist_out/`)
With a night, always `<name>_plan.csv` (status, target and standard times, airmasses, standard distance, Moon separation) and `<name>_visibility_N.png` (ten targets per plot).
Plus the formats chosen with `--formats`, each listing every target followed by its standard.
The default, `auto`, writes the chosen observatory's own format (IRTF → `irtf`, Keck → `keck`, Lick → `lick`, Paranal and La Silla → `eso`, Subaru → `subaru`, MMT → `mmt`, Las Campanas → `magellan`, La Palma → `cat`), and CSV otherwise:

| key | file | for |
| --- | --- | --- |
| `irtf` | `<name>_irtf.txt`, `<name>_irtf_nostds.txt` | IRTF starcat |
| `cat` | `cats/<name>_cat_N.cat` | ING object visibility (targets only, ten per file) |
| `txt` | `<name>_list.txt` | plain text table |
| `csv` | `<name>_list.csv` | spreadsheet / any software |
| `votable` | `<name>_list.vot` | TOPCAT |
| `keck` | `<name>_keck.txt` | Keck starlist (names up to 15 characters) |
| `lick` | `<name>_lick.txt` | Lick (Shane) starlist |
| `eso` | `<name>_eso_p2.txt` | ESO p2 *Import Target List* |
| `subaru` | `<name>_subaru_ope.txt` | Subaru OPE `<PARAMETER_LIST>` target lines |
| `mmt` | `<name>_mmt.cat` | MMT catalogue |
| `magellan` | `<name>_magellan.cat` | Magellan catalogue (rotator set to `0.0 OFF` and no guide stars: edit for your instrument) |

Each writer follows its observatory's published format.
Proper motions: the Keck, MMT and Magellan writers convert them to those telescopes' units; the Lick and ESO lists leave them out because their documentation does not say whether RA proper motion includes cos(Dec); the IRTF list copies the table's `pmra`/`pmdec` values unchanged, as version 1 did, because the units the IRTF StarCat expects could not be confirmed here, so check them against the StarCat user guide before relying on proper motions there.
Characters a format does not allow are removed from names (e.g. `* 42 Leo` becomes `42 Leo` for ESO).
Palomar and Gemini OT lists are not offered yet, as their target-list formats could not be confirmed from public documentation.

### Finder charts
Charts are cut from the CDS [hips2fits](https://alasky.cds.unistra.fr/hips-image-services/hips2fits) service (default survey 2MASS K, 2 arcmin), so no Java is needed.
Add another hips2fits service with `--hips2fits-url`.
If no hips2fits service answers and [Aladin Desktop](https://aladin.cds.unistra.fr/) is installed as an `aladin` command, the remaining charts are drawn in one Aladin session from the IRSA copies of the 2MASS, DSS2 and AllWISE surveys instead.

### Bundled data
`src/targetlist_generator/data/observatories.csv` is a snapshot of astropy's site registry (astropy-data `coordinates/sites.json`, BSD 3-Clause licence).

## Developing
```bash
pip install -e ".[dev,gui-test]"
pytest -m "not network"        # QT_QPA_PLATFORM=offscreen for the GUI tests on a headless machine
ruff check src tests
```
Version 1 (the single `targetlist_gen.py` script driving Aladin) is tagged `v1.0.4`.

## Citing
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6624861.svg)](https://doi.org/10.5281/zenodo.6624861)

The DOI above always points to the latest version; each release also has its own version DOI on Zenodo, which is best for citing the exact version you used.
Citation details are in `CITATION.cff`.
