"""Command line: ``targetlist make``, ``targetlist finder``, ``targetlist sites`` and ``targetlist gui``."""
import argparse
import sys

from . import __version__, finder, formats, io, sky, standards
from .pipeline import make_finders, make_targetlist

EPILOG = """examples:
  targetlist make targets.csv --date 2026-11-30 --site irtf
  targetlist make targets.fits --date 2026-11-30 --lat 28.76 --lon -17.88 --target-time-col exptime --finders
  targetlist make targets.csv --date 2026-11-30 --site lapalma --formats irtf,keck,votable --min-moon-sep off
  targetlist make targets.csv -r 20h 4h -d -30 60 --no-standards
  targetlist finder targets.csv J1023+2215 J0004-4044 --survey CDS/P/2MASS/J
  targetlist sites keck
  targetlist gui

Negative values work as plain numbers (-d -30 60); for negative sexagesimal strings use degrees instead."""


def formats_arg(value: str) -> list[str] | None:
    if value == 'auto':
        return None
    keys = [k.strip() for k in value.split(',') if k.strip()]
    unknown = sorted(set(keys) - set(formats.FORMATS))
    if unknown:
        raise argparse.ArgumentTypeError(f'unknown formats {unknown}; choose from {", ".join(formats.FORMATS)}')
    return keys


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog='targetlist', description=__doc__, epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    sub = parser.add_subparsers(dest='command', required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('target_list', help='CSV, FITS, ECSV, VOTable or IRTF .txt list; needs columns '
                                            'shortname (or name), ra, dec, spt and the magnitude column')
    common.add_argument('-o', '--out-dir', default='targetlist_out', help='output folder (default: %(default)s)')
    common.add_argument('-m', '--mag-col', default='tmassk', help='magnitude column (default: %(default)s)')
    common.add_argument('--survey', default=finder.SURVEYS[0], help='HiPS survey for finder charts, e.g. '
                        + ', '.join(finder.SURVEYS[:3]) + ' (default: %(default)s)')
    common.add_argument('--hips2fits-url', action='append', default=[],
                        help='extra hips2fits service tried before the CDS ones (repeatable), e.g. a mirror')
    common.add_argument('--fov', type=float, default=2.0, help='finder chart width [arcmin] (default: %(default)s)')

    make = sub.add_parser('make', parents=[common], help='cut and plan the list, write telescope target lists')
    night = make.add_argument_group('night (needs --date and --site or --lat/--lon)')
    night.add_argument('--date', help='evening date the night starts, YYYY-MM-DD')
    night.add_argument('--site', help="observatory name or alias, e.g. irtf, lapalma, keck (see 'targetlist sites')")
    night.add_argument('--lat', type=float, help='observatory latitude [deg], north positive')
    night.add_argument('--lon', type=float, help='observatory longitude [deg], east positive')
    night.add_argument('--twilight', type=float, default=-12.0, help='Sun altitude ending twilight [deg] '
                       '(default: %(default)s; -18 for astronomical)')
    night.add_argument('--max-airmass', type=float, default=2.0, help='(default: %(default)s)')
    night.add_argument('--min-time', type=float, default=30.0,
                       help='minutes a target must spend below --max-airmass (default: %(default)s)')
    night.add_argument('--min-moon-sep', default='auto',
                       help="minimum Moon separation: auto (10 + 30 x illumination deg), off, or degrees "
                            "(default: %(default)s)")
    times = make.add_argument_group('time per block (target, calibrations, standard) [min]')
    times.add_argument('--target-time', type=float, default=30.0, help='same for every target (default: %(default)s)')
    times.add_argument('--target-time-col', help='table column with each target\'s minutes instead')
    times.add_argument('--cal-time', type=float, default=5.0, help='calibrations per target (default: %(default)s)')
    times.add_argument('--std-time', type=float, default=10.0, help='standard per target (default: %(default)s)')
    box = make.add_argument_group('RA/Dec box (optional, combines with the night cut)')
    box.add_argument('-r', '--ra-range', nargs=2, metavar=('LOWER', 'UPPER'),
                     help='degrees or hours like 20h; wraps through 0 when LOWER > UPPER')
    box.add_argument('-d', '--dec-range', nargs=2, metavar=('LOWER', 'UPPER'), help='degrees or e.g. 15d30m')
    std = make.add_argument_group('telluric standards')
    std.add_argument('--no-standards', action='store_true', help='skip pairing standards')
    std.add_argument('--standards-file', help='CSV like the bundled list (shortname, ra, dec, spt, band columns)')
    std.add_argument('--std-band', default='K', choices=standards.BANDS, help='(default: %(default)s)')
    std.add_argument('--std-spt', help='only standards of this type, e.g. A0V (default: B9V-A1V)')
    std.add_argument('--std-airmass-tol', type=float, default=0.1,
                     help='standard within this airmass of the target when possible (default: %(default)s)')
    std.add_argument('--std-max-sep', type=float, default=20.0,
                     help='standard within this many deg of the target when possible (default: %(default)s)')
    make.add_argument('--formats', type=formats_arg, default=None,
                      help=f'comma-separated output formats from {", ".join(formats.FORMATS)} '
                           '(default: auto = the observatory\'s own format, else csv)')
    make.add_argument('--finders', action='store_true', help='also make finder charts for every listed target')

    fc = sub.add_parser('finder', parents=[common], help='finder charts for named targets in the list')
    fc.add_argument('names', nargs='+', help="shortnames, or 'all'")

    sites = sub.add_parser('sites', help='list the bundled observatories (optionally matching a search)')
    sites.add_argument('search', nargs='?', default='', help='part of a name or alias')

    sub.add_parser('gui', help='open the graphical interface')
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == 'gui':
        try:
            from .gui import main as gui_main
        except ImportError as err:  # PySide6 is an optional extra
            sys.exit(f'targetlist: the GUI needs PySide6 ({err}); install "targetlist-generator[gui]"')
        gui_main()
        return
    if args.command == 'sites':
        obs = sky.observatories()
        hit = (obs['name'] + ';' + obs['aliases']).str.contains(args.search, case=False, regex=False)
        print(obs[hit].to_string(index=False))
        return
    urls = tuple(args.hips2fits_url) + finder.HIPS2FITS
    try:
        if args.command == 'make':
            make_targetlist(args.target_list, out_dir=args.out_dir, mag_col=args.mag_col,
                            target_time=args.target_time, target_time_col=args.target_time_col,
                            cal_time=args.cal_time, with_standards=not args.no_standards,
                            standards_file=args.standards_file, std_band=args.std_band, std_spt=args.std_spt,
                            std_time=args.std_time, std_airmass_tol=args.std_airmass_tol,
                            std_max_sep=args.std_max_sep, outputs=args.formats, finders=args.finders,
                            survey=args.survey, fov=args.fov, urls=urls, date=args.date, site=args.site,
                            lat=args.lat, lon=args.lon, twilight=args.twilight, max_airmass=args.max_airmass,
                            min_time=args.min_time, min_moon_sep=args.min_moon_sep,
                            ra_range=args.ra_range, dec_range=args.dec_range)
        else:
            targets = io.read_targets(args.target_list, args.mag_col)
            names = targets['shortname'] if args.names == ['all'] else args.names
            make_finders(targets, names, args.out_dir, args.survey, args.fov, urls)
    except (ValueError, FileNotFoundError, RuntimeError) as err:  # user-facing problems: message, no traceback
        sys.exit(f'targetlist: error: {err}')
