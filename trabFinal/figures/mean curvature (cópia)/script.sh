# FILES="./"
for f in *.png
do
	echo "Processiong $f"
	convert -crop 1200x400+0+250 "$f" "$f"
done