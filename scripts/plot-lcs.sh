starinfo=""
energyout=""
arrivalout=""
flareout=""
starperiod=2.69

files = `ls ./*.fits`
for f in $files
do
    lightcurve_analysis.py f $starperiod > $flareout
    #correct_flares.py flareout > flareout
    flare_energy.py $flareout $starinfo >> $energyout
    flare_interarrival.py $flareout >> $arrivalout
done
