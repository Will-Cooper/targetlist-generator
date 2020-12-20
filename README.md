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
