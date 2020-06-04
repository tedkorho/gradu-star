for f in ./*/*/*
do
    dir="${f%/*}"
    python3 trend_curve.py "$f" 2.59 --plot > "$dir/trended_lightcurve.out"
    echo "analyzing $dir"
done
