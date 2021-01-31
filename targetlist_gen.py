"""
A script to take an initial target list and cut it to a given visibility range, find standards
and generate finder charts.

You should provide EITHER -t and -r and -d AND/OR -n as optional parameters.
The former generates the target list whilst the latter provides a finder chart.

If you are giving a negative value in a string (e.g. -15h), argparse has a bug which will fail the script.
You will have to give a literal string on the command line (e.g. " -15h").
Beware if not giving an RA limit and relying on the conversion from LST, the sunset and sunrise times you
give MUST be UTC not local time.
Input files which are not .csv or .fits will probably fail unless one has provided some exact columns.

Methods
-------
live_print
    Printing to stdout in a live fashion
df_editor
    Edit dataframe column names and positional data types
ang_converter
    Converts angles from degrees to sexadecimal
tab_parser
    Checks the type of input table, e.g. fits, csv, txt
find_ra_lims
    Determines the RA limits based on UTC sunset/ sunrise and observatory co-ordinates
find_dec_lims
    Determines the DEC limits based on latitude +/- 30 degrees
main
    Main module, handles arguments and potential errors

Classes
-------
TableEditor
    Creating cut target list in both .txt form and IRTF form
FinderCharts
    Construct finder chart for a given object
"""
# standard library imports
import argparse
import glob
import os
import sys
from typing import Union, Tuple, List
import warnings

# external libraries
from astropy.coordinates import Angle, SkyCoord, EarthLocation, Latitude
import astropy.units as u
from astropy.table import Table
from astropy.time import Time
import numpy as np
import pandas as pd


def live_print(s: str, end: str = '\n'):
    """
    Writes to stdout in live fashion unlike inbuilt print

    Parameters
    ----------
    s
        String to be written
    end
        Line end of string
    """
    sys.stdout.write(s + end)  # write to stdout
    sys.stdout.flush()  # flush stdout pending outputs
    return


def df_editor(df: pd.DataFrame) -> pd.DataFrame:
    """
    Edits dataframe column names and datatype of positional columns

    Parameters
    ----------
    df
        Dataframe of objects

    Returns
    -------
        Dataframe of edited objects
    """
    cols = df.columns  # columns of dataframe
    df = df.rename(columns={i: i.lower() for i in cols})  # port all columns into lower case
    df = df.rename(columns={'name': 'shortname', 'median_sptname': 'spt',
                            'gaiag_spt': 'spt'})  # edit these specific columns if they exist
    if df['ra'].dtype != 'float64':  # if RA column is sexadecimal
        ra = df['ra'].values
        if ':' in ra[0]:  # if ra in form hh:mm:ss.ss
            for j, i in enumerate(ra):
                i = i.split(':')
                ra[j] = "".join([i[0], 'h', i[1], 'm', i[2], 's'])  # convert to hms
        ra = Angle(ra)  # create as angles
        df['ra'] = ra.deg  # convert to degrees
    if df['dec'].dtype != 'float64':  # if Dec column is sexadecimal
        dec = df['dec'].values
        if ':' in dec[0]:  # if dec in form dd:mm:ss.s
            for j, i in enumerate(dec):
                i = i.split(':')
                dec[j] = "".join([i[0], 'd', i[1], 'm', i[2], 's'])  # convert to dms
        dec = Angle(dec)  # create as angles
        df['dec'] = dec.deg  # convert to degrees
    return df


def ang_converter(ra: Angle, dec: Angle) -> Tuple[int, int, float, str, int, float]:
    """
    Convert ra and dec from degrees to sexadecimal

    Parameters
    ----------
    ra
        RA of object in degrees
    dec
        Dec of object in degrees

    Returns
    -------
    rah
        Integer hour of RA
    ram
        Integer minutes of RA
    ras
        Float seconds of RA, rounded to 2dp
    decd
        Integer degrees of Dec (including sign)
    decm
        Integer minutes of Dec
    decs
        Float seconds of Dec, rounded to 1dp
    """
    rah = ra.hms[0].astype(int)  # ra hours
    ram = ra.hms[1].astype(int)  # ra mins
    ras = round(ra.hms[2], 2)  # ra seconds
    decd = (dec.signed_dms[0] * dec.signed_dms[1]).astype(int)  # dec degrees
    decm = dec.signed_dms[2].astype(int)  # dec mins
    decs = round(dec.signed_dms[3], 1)  # dec seconds
    return rah, ram, ras, decd, decm, decs


def spt_from_num(sptnum: Union[float, int]) -> str:
    """
    Converts spectral type number to spectral type

    Parameters
    ----------
    sptnum
        Spectral type number of object

    Returns
    -------
        Spectral type of object
    """
    def spt_make(sptype: str, sub: int) -> str:
        """
        Joins strings of spectral type and sub spectral type

        Parameters
        ----------
        sptype
            Main spectral type of object
        sub
            Sub type
        Returns
        -------
        _
            Joined up string of main type and sub type
        """
        def round_val() -> str:
            """
            Rounds values to nearest 0.5

            Returns
            -------
            _
                String of rounded value
            """
            val = sptnum - sub  # difference between spectral type number of object and presumed sub type
            if val >= 10:  # if that difference is >= 10, spectral type number is earlier than M or later than Y
                raise ArithmeticError('Spectral type number incorrect')
            rval = 0.5 * round(val / 0.5)  # round to nearest 0.5
            return str(rval)
        return "".join([sptype, round_val()])
    if sptnum < 70:
        spt = spt_make('M', 60)  # M dwarf
    elif sptnum < 80:
        spt = spt_make('L', 70)  # L dwarf
    elif sptnum < 90:
        spt = spt_make('T', 80)  # T dwarf
    else:
        spt = spt_make('Y', 90)  # Y dwarf
    if '.0' in spt:  # prettify interger spectral types
        return spt[:2]
    else:
        return spt


def tab_parser(tab: str):
    """
    Handles different file types differently, tries csv, then fits, then fallback is txt

    Parameters
    ----------
    tab
        Name of table to be opened into pandas dataframe
    Returns
    -------
    df
        Pandas dataframe post editing
    """
    if 'csv' in tab:
        df = df_editor(pd.read_table(tab, sep=','))  # load and edit table
    elif 'fits' in tab:
        df = df_editor(Table.read(tab).to_pandas())
        df['shortname'] = [i.decode('utf-8').strip() for i in df['shortname']]
    else:
        df = df_editor(pd.read_table(tab, sep='\t', names=['index', 'name',
                                                           'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', 'tmassk'],
                                     usecols=list(range(0, 9))))  # load and edit table if using .txt file
    return df


def find_ra_lims(sunset: str, sunrise: str, longitude: List[int], latitude: List[int]):
    """
    Determines the RA limits based on the local sidereal time

    Parameters
    ----------
    sunset
        UTC sunset time
    sunrise
        UTC sunrise time
    longitude
        Longtiude of observatory
    latitude
        Latitude of observatory

    Returns
    -------
    ramin
        RA start time
    ramax
        RA end time
    """
    obs_loc = EarthLocation(lat=(*latitude, ), lon=(*longitude, ))
    tset = Time(sunset, scale='utc', location=obs_loc)
    trise = Time(sunrise, scale='utc', location=obs_loc)
    ramin = str(np.floor(tset.sidereal_time('apparent').value)) + 'h'
    ramax = str(np.ceil(trise.sidereal_time('apparent').value)) + 'h'
    return ramin, ramax


def find_dec_lims(latitude: List[int]):
    """
    Grabs approximate declination limits based on the latitude +/- 30 degrees

    Parameters
    ----------
    latitude
        Latitude of the observatory

    Returns
    -------
    decmin
        Lower declination limit
    decmax
        Upper declination limit
    """
    lat = Latitude(f'{latitude[0]}d{latitude[1]}m')
    decmin = np.floor(lat - 30 * u.deg).value
    decmax = np.ceil(lat + 30 * u.deg).value
    if decmin < -90:
        decmin = 90
    if decmax > 90:
        decmax = 90
    return decmin, decmax


class TableEditor:
    """
    Creates cut target list in useful formats

    Methods
    -------
    limits_transform
        Converts limits from sexadecimal to degrees
    df_cut
        Cut dataframe on ra & dec limits
    standards_load
        Load the standards fits file and edit to useful form
    closest_standard
        Find the closest standard for each object in target list
    df_connect
        Create dataframe of alternating object and respective standard
    fwrite
        Write file in IRTF format
    """

    def __init__(self, tab: str, std_tab: str, *limits):
        """
        Constructor method for class, all functions handled from here

        Parameters
        ----------
        tab
            Name of target list file
        std_tab
            Name of standard list file
        limits
            RA & Dec lower & upper limits
        """
        df = tab_parser(tab)
        live_print(f'{tab} loaded')
        self.ramin, self.ramax, self.decmin, self.decmax = self.limits_transform(limits)  # ensure limits are in degs
        live_print(f'Limits understood -- RA: {self.ramin:.2f} - {self.ramax:.2f} degrees'
                   f' and Dec: {self.decmin:.1f} - {self.decmax:.1f} degrees')
        df['pmra'] = [0.0, ] * len(df)  # set pmra to 0
        df['pmdec'] = [0.0, ] * len(df)  # set pmdec to 0
        self.df = self.df_cut(df)  # cut on ra & dec limits
        self.df.reset_index(drop=True, inplace=True)
        live_print(f'There are {len(self.df)} targets')
        self.dfstd = self.df_cut(self.standards_load(std_tab))  # load the standards
        live_print('Cut on limits and standards loaded')
        self.dfc = self.closest_standard()  # find the closest standards
        live_print('Closest standards found')
        self.dfcomb = self.df_connect()  # create dataframe of alternating object & respective standard
        live_print('Tables concatenated')
        self.fnameout = tab[:tab.find('.')] + '_edited.txt'  # name of file to be written
        self.dfcomb.to_csv(self.fnameout, header=False, sep='\t')  # write file
        self.fwriter()  # write file in IRTF form
        self.cat_writer()  # write file in ING form
        live_print('File written')
        return

    @staticmethod
    def limits_transform(*limits) -> np.ndarray:
        """
        Transfer limits from sexadecimal to degrees if not already done so

        Parameters
        ----------
        limits[0]
            Tuple of limits: ra lower, ra upper, dec lower, dec upper

        Returns
        -------
        limout
            Tuple of limits: ra lower, ra upper, dec lower, dec upper; all in degrees
        """
        limout = np.empty(4)  # empty array
        for i, var in enumerate(limits[0]):  # over every limit
            if type(var) == str:  # check if it is in sexadecimal form
                var = Angle(var).deg  # convert to degrees
            limout[i] = var  # assign to new array
        return limout

    def df_cut(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cut the target list on ra & dec limits
        Parameters
        ----------
        df
            Dataframe of targetlist

        Returns
        -------
        df
            Dataframe of targetlist, now cut
        """
        if self.ramax > self.ramin:
            df = df.loc[np.logical_and(df['ra'] > self.ramin, df['ra'] < self.ramax),  # cut on ra
                        ['shortname', 'ra', 'dec', 'pmra', 'pmdec', 'spt', 'tmassk']]  # only keep these columns
        else:
            df = df.loc[np.logical_or(df['ra'] > self.ramin, df['ra'] < self.ramax),  # cut on ra
                        ['shortname', 'ra', 'dec', 'pmra', 'pmdec', 'spt', 'tmassk']]  # only keep these columns
        df = df.loc[np.logical_and(df['dec'] > self.decmin, df['dec'] < self.decmax)]  # cut on dec
        df['epoch'] = [2000, ] * len(df)
        if df['spt'].dtype in ('float64', 'int64'):
            df['spt'] = [spt_from_num(i) for i in df['spt'].values]  # convert spectral type numbers to spectral type
        if len(df) < 1:  # check we haven't cut to an empty target list
            raise ArithmeticError('Check limits, cutting has created an empty table')
        df = df[['shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', 'tmassk']]  # set column order
        return df

    @staticmethod
    def standards_load(f: str) -> pd.DataFrame:
        """
        Load the standards file and convert into a useful form

        Parameters
        ----------
        f
            The filename of the standards file

        Returns
        -------
        dfstd
            Dataframe of standards
        """
        dfstd = pd.read_table(f, sep=',')
        dfstd = dfstd.rename(columns={'K': 'tmassk'})  # rename columns to be same as target dataframe
        dfstd = dfstd[['shortname', 'ra', 'dec', 'spt', 'tmassk']]  # only take these columns from dataframe
        dfstd['pmra'] = [0.0, ] * len(dfstd)
        dfstd['pmdec'] = [0.0, ] * len(dfstd)
        dfstd['shortname'] = [i.strip().replace(' ', '_') for i in dfstd['shortname'].values]
        dfstd = dfstd.iloc[np.flatnonzero(['HD' in std for std in dfstd['shortname']])]
        return dfstd

    def closest_standard(self) -> pd.DataFrame:
        """
        Determines the closest standard to each object

        Returns
        -------
        dfc
            Dataframe of the closest standards
        """
        dfc = pd.DataFrame(columns=['shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', 'tmassk'])  # dataframe
        live_print('Finding closest standards', end='')
        rastds = self.dfstd['ra'].values * u.deg
        decstds = self.dfstd['dec'].values * u.deg
        cstds = SkyCoord(rastds, decstds)  # coordinate of standards
        for i in range(len(self.df)):  # over every target row
            ra = self.df.iloc[i]['ra'] * u.deg + 11.25 * u.deg  # RA of target + 45m
            dec = self.df.iloc[i]['dec'] * u.deg  # Dec of target
            cobj = SkyCoord(ra, dec)  # coordinate of object
            dist = cobj.separation(cstds)  # distance between object and standards
            ind = dist.argmin()  # index of minimum difference
            row = self.dfstd.iloc[[ind]]  # the row in the standards dataframe that is closest
            dfc = pd.concat([dfc, row])  # append that row to new dataframe
            live_print('.', '')
        live_print('')
        return dfc

    @staticmethod
    def coord_conv(ra: float, dec: float):
        """
        Constructs string of ra and dec in sexadecimal from degrees

        Parameters
        ----------
        ra
            RA of target in degrees
        dec
            Dec of target in degrees
        Returns
        -------
        rastring
            RA of target in form hh:mm:ss.ss
        decstring
            Dec of target in form dd:mm:ss.s
        """
        ra = Angle(ra * u.deg)
        dec = Angle(dec * u.deg)
        rah, ram, ras, decd, decm, decs = ang_converter(ra, dec)  # convert and unpack ra & dec
        rastring = f'{rah}:{ram}:{ras}'
        decstring = f'{decd}:{decm}:{decs}'
        return rastring, decstring

    def df_connect(self) -> pd.DataFrame:
        """
        Connects targets with their closest standards in one dataframe

        Returns
        -------
        dfcomb
            Combined dataframe of targets and respective standards
        """
        def row_editor(row: np.ndarray) -> np.ndarray:
            """
            Edit each row by converting and rounding

            Parameters
            ----------
            row
                Row of data to be edited
            Returns
            -------
            row
                Edited row
            """
            row[1], row[2] = self.coord_conv(row[1], row[2])  # convert ra & dec to hms/ dms
            row[-1] = round(row[-1], 2)  # round the mag column
            row[4] = round(row[4], 2)  # round the pmra
            row[5] = round(row[5], 2)  # round the pmdec
            return row

        dfcomb = pd.DataFrame(columns=['shortname', 'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', 'tmassk'])
        live_print('Combining tables', '')
        for i in range(len(self.df)):  # over every row of dataframe
            rowobj = row_editor(np.array(self.df.iloc[[i]])[0])  # gather row of object and edit
            rowstd = row_editor(np.array(self.dfc.iloc[[i]])[0])  # gather row of standard and edit
            # append first the object then the standard as rows to the same table
            rowobj = pd.DataFrame(rowobj.reshape(-1, len(rowobj)))
            rowstd = pd.DataFrame(rowstd.reshape(-1, len(rowstd)))
            dfcomb = pd.concat([dfcomb, rowobj, rowstd], ignore_index=True)
            live_print('.', end='')
        live_print('')
        dfcomb.index += 1  # dataframes are indexed from 0, we need to start from 1
        return dfcomb

    def fwriter(self):
        """
        Writes file in form understandable by starcat/ IRTF
        """
        with open(self.fnameout[:self.fnameout.find('.')] + '_irtf', 'w+') as f:  # open irtf file as write
            # header for IRTF
            f.write('#+-----+-----+-----+-----+-----+-----+-----+-----+-------+\n')
            f.write('#|Index| Name|  RA | DEC |Epoch| PMRA|PMDEC|SpT  | K2MASS|\n')
            f.write('#+-----+-----+-----+-----+-----+-----+-----+-----+-------+\n')
            with open(self.fnameout[:self.fnameout.find('.')] + 'tmp', 'w') as ftmp:
                with open(self.fnameout, 'r') as fin:  # open normal txt target list
                    for line in fin:  # append every line to the irtf file
                        line = line.rstrip('\n').rstrip('\t') + '\n'  # fixing a write issue
                        ftmp.write(line)
                        f.write(line)
        os.rename(self.fnameout[:self.fnameout.find('.')] + 'tmp', self.fnameout)
        with open(self.fnameout[:self.fnameout.find('.')] + '_irtf', 'r+') as fin:  # open irtf file as read
            with open(self.fnameout[:self.fnameout.find('.')] + '_irtf_nostds', 'w+') as fout:  # no stds file
                for line in fin:
                    if np.all([spt not in line for spt in ('B9V', 'A0V', 'A1V')]):
                        fout.write(line)
        return

    def cat_writer(self):
        """
        Writes in the .cat format (for ING object visibility)
        """
        j, length = 0, len(self.df)
        chunk = 10
        if chunk > len(self.df):
            chunk = len(self.df)
        while j < length:
            top = j + chunk
            if top > length:
                top = length
            df = self.df.loc[j: top, 'shortname': 'epoch'].copy()
            if len(df) == 0:
                break
            ra, dec = [], []
            for i in range(len(df)):
                converted = self.coord_conv(np.array(df.iloc[[i]]['ra'])[0], np.array(df.iloc[[i]]['dec'])[0])
                ra.append(converted[0])
                dec.append(converted[1])
            df['ra'] = ra
            df['dec'] = dec
            df.to_csv('cats/' + self.fnameout[:self.fnameout.find('.')] + f'_cat_{j // 10 + 1}.cat',
                      header=False, sep=' ', index=False)
            j += chunk
        return


class FinderCharts:
    """
    Constructs finder chart for given object and band

    Methods
    -------
    find_coords
        Finds the coordinates in degrees from the given table
    query_aladin
        Constructs and sends query to Aladin
    """

    def __init__(self, tab: str, name: str, fov: Union[float, int],
                 band: str, aladin: str, view: bool):
        """
        Constructor method of class, all functions handled from here

        Parameters
        ----------
        tab
            Name of target list file to look in
        name
            Name of target to look for
        fov
            Field of view to make finder chart in
        band
            Survey/ band that Aladin can understand
        aladin
            User path to Aladin.jar file
        view
            Optional argument for viewing data
        """
        if 'csv' in tab:
            df = df_editor(pd.read_table(tab, sep=','))  # load and edit table
        else:
            df = df_editor(pd.read_table(tab, sep='\t', names=['index', 'name',
                                                               'ra', 'dec', 'epoch', 'pmra', 'pmdec', 'spt', 'tmassk'],
                                         usecols=list(range(0, 9))))  # load and edit table if using .txt file
        self.name = name
        if self.name not in df['shortname'].values:  # check target is in target list table
            raise ValueError(f'{self.name} not in {tab}')
        self.df = df
        self.ra, self.dec = self.find_coords()  # get the coordinates in degrees
        self.fov = fov
        self.band = band
        self.aladin_path = aladin
        self.view = view
        self.query_aladin()  # query Aladin
        return

    def find_coords(self) -> Tuple[str, str]:
        """
        Construct strings from ra & dec which are in degrees

        Returns
        -------
        rastring
            String of ra position in form hh mm ss.ss
        decstring
            String of dec position in form dd mm ss.s
        """
        ra = Angle(self.df.loc[self.df['shortname'] == self.name, 'ra'].values[0] * u.deg)  # get ra from table
        dec = Angle(self.df.loc[self.df['shortname'] == self.name, 'dec'].values[0] * u.deg)  # get dec from table
        rah, ram, ras, decd, decm, decs = ang_converter(ra, dec)  # convert degrees to hms & unpack
        rastring = f'{rah} {ram} {ras}'  # construct string
        decstring = f'{decd} {decm} {decs}'  # construct string
        return rastring, decstring

    def query_aladin(self):
        """
        Constructs query and transports to aladin for creation of finder chart
        """
        pos_qry = f'get hips({self.band}) {self.ra} {self.dec} {self.fov}\''  # query survey
        sav = f'save 524x460 fcharts/{self.name}_{self.fov}.png'  # save function
        aladin = f'java -jar {self.aladin_path}'  # get correct aladin path
        os.system("echo \"" + pos_qry + "; sync; " + sav + "; quit\" | " + aladin + " -nogui")  # construct
        if self.view:  # if user wants to view finder chart
            os.system("echo \"" + pos_qry + "; sync\" | " + aladin)  # construct for viewing
            live_print('')
        return


def main():
    """
    Main module handling argument parsing and error handling

    """
    def floatorstr(var: str) -> Union[float, str]:
        """
        Determine whether a variable is either a float or string from a string

        Parameters
        ----------
        var
            Variable is either float or string but curremt;y string

        Returns
        -------
        var
            Either float or string variable
        """
        try:
            return float(var)
        except ValueError:
            return var

    # argument parsing
    myargs = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    myargs.add_argument('-t', '--target-list', help='Name of target list', required=True)
    tgtlist = myargs.add_argument_group('Target List')
    tgtlist.add_argument('-c', '--create-list', help='Select if you want to make the target list',
                         action='store_true', default=False)
    tgtlist.add_argument('-r', '--ra-range', nargs=2, metavar=('RA_LOWER_LIMIT', 'RA_UPPER_LIMIT'),
                         help='RA range in form "000.0" or "00h00m00.0s"', default=(0, 360))
    tgtlist.add_argument('-d', '--dec-range', nargs=2, metavar=('DEC_LOWER_LIMIT', 'DEC_UPPER_LIMIT'),
                         help='Dec range in form "00.0" or "s00d00m00s"', default=(-90, 90))
    tgtlist.add_argument('-s', '--standard_list', help='Name of standard list',
                         default='large_standards_bright.csv')
    fchart = myargs.add_argument_group('Finder Chart')
    fchart.add_argument('-n', '--name', help='Name of target')
    fchart.add_argument('-b', '--aladin-band', default='CDS/P/2MASS/K', help='Band to get finder chart for')
    fchart.add_argument('-f', '--fov', type=int, default=1, help='Field of View in arcminutes for finder chart')
    fchart.add_argument('-a', '--aladin-path', help='Full path to aladin jar file', default='~/.Aladin/Aladin.jar')
    fchart.add_argument('-v', '--view-target', action='store_true', default=False,
                        help='Select to view object in Aladin')
    findlims = myargs.add_argument_group('Find RA Limits')
    findlims.add_argument('--find-lims', help='Find RA limits from Sunset/ Sunrise', default=False, action='store_true')
    findlims.add_argument('--sunset', help='UTC Sunset Time in form YYYY-MM-DDTHH:MM (T as written)')
    findlims.add_argument('--sunrise', help='UTC Sunrise Time in form YYYY-MM-DDTHH:MM (T as written)')
    findlims.add_argument('--longitude', help='Observatory longitude in form DDD MM, e.g. -153 28', type=int, nargs=2,
                          metavar=('LONG_DEG', 'LONG_MIN'))
    findlims.add_argument('--latitude', help='Observatory latitude in form DD MM, e.g. 19 49', type=int, nargs=2,
                          metavar=('LAT_DEG', 'LAT_MIN'))
    args = myargs.parse_args()
    do_list = bool(args.create_list)
    do_view = bool(args.view_target)
    tab = str(args.target_list)
    tgt = args.name
    rarange = tuple(args.ra_range)
    decrange = tuple(args.dec_range)
    std = args.standard_list
    band = str(args.aladin_band)
    fov = args.fov
    aladin_path = args.aladin_path
    ramin, ramax = floatorstr(rarange[0]), floatorstr(rarange[1])
    decmin, decmax = floatorstr(decrange[0]), floatorstr(decrange[1])
    do_lims = bool(args.find_lims)
    sunset = args.sunset
    sunrise = args.sunrise
    longitude = args.longitude
    latitude = args.latitude

    # error handling
    if len(glob.glob('fcharts/')) < 1 or len(glob.glob('cats/')) < 1:  # if finder chart folder missing
        class MakingDirWarning(UserWarning):
            pass
        warnings.warn("Attempting to create folders 'fcharts' and 'cats' (will fail later if no permissions)",
                      MakingDirWarning, stacklevel=3)
        os.mkdir('fcharts')  # attempt to make it
        os.mkdir('cats')
    if len(glob.glob(std)) < 1 and tgt is None:  # look for standards file
        raise FileNotFoundError(f'Require standards file {std} when generating target list')
    if len(glob.glob(tab)) < 1:  # look for stated table
        raise FileNotFoundError(f'Cannot find {tab}')
    else:  # otherwise proceed
        df = tab_parser(tab)
        cols = list(df.columns)
        if np.any([col not in cols for col in ('ra', 'dec', 'shortname', 'tmassk', 'spt')]):  # need specifc columns
            raise ValueError(f'Need columns "ra", "dec", "tmassk", "spt"'
                             f' and "name" or "shortname" in table {tab} (any case)')

    # running selected tasks
    if tgt is not None:  # if finder chart requested
        tgt = str(tgt)
        FinderCharts(tab, name=tgt, fov=fov, band=band, aladin=aladin_path, view=do_view)  # finder chart for object
    if do_list:  # if creation of list requested
        # some warnings if using default ra & dec ranges
        class DefaultVarsWarning(UserWarning):
            pass
        if 0 == ramin and 360 == ramax and not do_lims:
            warnings.warn("Using default values of RA range (0 - 360)", DefaultVarsWarning, stacklevel=3)
        if -90 == decmin and 90 == decmax and latitude is None:
            warnings.warn("Using default values of Dec range (-90 - 90)", DefaultVarsWarning, stacklevel=3)
        elif -90 == decmin and 90 == decmax:
            print('Approximating DEC limits based on observatory latitude +/- 30 degrees')
            decmin, decmax = find_dec_lims(latitude)
        if 0 == ramin and 360 == ramax and do_lims:
            if np.any([val is None for val in (sunrise, sunset, longitude, latitude)]):
                print("Need to provide sunset, sunrise, longitude and latitude to find RA limits")
                warnings.warn("Falling back to default values of RA range (0 - 360)", DefaultVarsWarning, stacklevel=3)
            else:
                print('Approximating RA limits based on night time and coordinates')
                ramin, ramax = find_ra_lims(sunset, sunrise, longitude, latitude)
        TableEditor(tab, std, ramin, ramax, decmin, decmax)  # create table using cuts from given limits
    elif not do_list and tgt is None:  # otherwise nothing is going to happen so warn about it
        class InvalidChoiceWarning(UserWarning):
            pass
        warnings.warn("Either select a target by name (-n) or create the target list (-c)",
                      InvalidChoiceWarning, stacklevel=3)
        myargs.parse_args(['-h', '-t'])
    return


if __name__ == '__main__':  # if running as script
    main()  # call main module
