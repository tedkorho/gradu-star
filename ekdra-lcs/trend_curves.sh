flareout="temp.out"
energyout="energies_ekdra.out"
arrivalout="interarrival_times_ekdra.out"
flareout2="temp2.out"
flareout3="temp3.out"

for f in ./*/*/*/*
do
    dir="${f%/*}"
    echo "analyzing $dir"
	python3 lightcurve_analysis.py "$f" 2.59 --plot > "$flareout"
	cp "$flareout" "$dir/lightcurve_data.out"
	python3 trend_iteration.py "$dir/lightcurve_data.out" 2.59 --plot > "$flareout2"
    python3 trend_iteration.py "$flareout2" 2.59 --plot > "$flareout3"
	cp "$flareout3" "$dir/lightcurve_data.2.out"
#	python3 flare_energy.py "$flareout2" "$f" >> "$energyout"    
#	python3 flare_interarrival.py "$flareout2" >> "$arrivalout"
    rm "$flareout"
	rm "$flareout2"
	rm "$flareout3"
done

# ADD --show to the trend_iteration.py commands to view results during calculation
