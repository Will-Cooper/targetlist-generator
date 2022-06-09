# Targetlist Generator

## Installation
Ensure dependencies are correct, see `requirements.txt`; best course of action is to use conda:
```bash
conda create --name targetlist --file requirements.txt
conda activate targetlist
```
Then clone or fork this repo:
```bash
git clone https://github.com/Will-Cooper/targetlist-generator.git
cd targetlist-generator/
```

## Running
A script to automatically cut target lists into RA and DEC limits, where RA can be calculated just off sunset/ sunrise times (in UTC!) if desired.  
Run in terminal with 
```bash
python targetlist_gen.py -h
``` 
to get the full set of command line arguments.  

Can handle .fits or .csv, maybe .txt if the format reads nicely.
Will return a cut list in the form expected at IRTF alongside multiple .cat files in cats/. This cut list will be in two forms, one with just the targets (...nostds) and the other interspersed with the closest standard star, assuming +30 mins from target zenith.  

One can also use the script to generate a finder chart via [aladin](https://aladin.u-strasbg.fr/ "Aladin") (note the option --aladin-path).  
At current there are some specific expected columns, most notably: 'ra', 'dec', 'shortname', 'tmassk', 'spt'. RA and DEC can be given either sexidecimally or decimalised. Can possibly be changed if different mags are desired but K is required for IRTF.

```text
usage: targetlist_gen.py [-h] -t TARGET_LIST [-c]
                         [-r RA_LOWER_LIMIT RA_UPPER_LIMIT]
                         [-d DEC_LOWER_LIMIT DEC_UPPER_LIMIT]
                         [-s STANDARD_LIST] [-n NAME] [-b ALADIN_BAND]
                         [-f FOV] [-a ALADIN_PATH] [-v] [--find-lims]
                         [--sunset SUNSET] [--sunrise SUNRISE]
                         [--longitude LONG_DEG LONG_MIN]
                         [--latitude LAT_DEG LAT_MIN]

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

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET_LIST, --target-list TARGET_LIST
                        Name of target list

Target List:
  -c, --create-list     Select if you want to make the target list
  -r RA_LOWER_LIMIT RA_UPPER_LIMIT, --ra-range RA_LOWER_LIMIT RA_UPPER_LIMIT
                        RA range in form "000.0" or "00h00m00.0s"
  -d DEC_LOWER_LIMIT DEC_UPPER_LIMIT, --dec-range DEC_LOWER_LIMIT DEC_UPPER_LIMIT
                        Dec range in form "00.0" or "s00d00m00s"
  -s STANDARD_LIST, --standard_list STANDARD_LIST
                        Name of standard list

Finder Chart:
  -n NAME, --name NAME  Name of target
  -b ALADIN_BAND, --aladin-band ALADIN_BAND
                        Band to get finder chart for
  -f FOV, --fov FOV     Field of View in arcminutes for finder chart
  -a ALADIN_PATH, --aladin-path ALADIN_PATH
                        Full path to aladin jar file
  -v, --view-target     Select to view object in Aladin

Find RA Limits:
  --find-lims           Find RA limits from Sunset/ Sunrise
  --sunset SUNSET       UTC Sunset Time in form YYYY-MM-DDTHH:MM (T as
                        written)
  --sunrise SUNRISE     UTC Sunrise Time in form YYYY-MM-DDTHH:MM (T as
                        written)
  --longitude LONG_DEG LONG_MIN
                        Observatory longitude in form DDD MM, e.g. -153 28
  --latitude LAT_DEG LAT_MIN
                        Observatory latitude in form DD MM, e.g. 19 49
  ```

### Citing
[![DOI](https://zenodo.org/badge/323075587.svg)](https://zenodo.org/badge/latestdoi/323075587)
