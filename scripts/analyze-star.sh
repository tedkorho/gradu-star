flareout="foo.dat"
starperiod=2.69

files = `ls ./*/*.fits`
for f in $files
do
    lightcurve_analysis.py f $starperiod > $flareout --plot
done
