#!/bin/bash
mylist="PDBDEV_00000001_resize, PDBDEV_00000002_resize, PDBDEV_00000003_resize, PDBDEV_00000004_resize, PDBDEV_00000005_resize, PDBDEV_00000006_resize, PDBDEV_00000007_resize, PDBDEV_00000008_resize, PDBDEV_00000009_resize, PDBDEV_00000010_resize, PDBDEV_00000011_resize, PDBDEV_00000012_resize, PDBDEV_00000014_resize, PDBDEV_00000015_resize, PDBDEV_00000016_resize, PDBDEV_00000017_resize, PDBDEV_00000018_resize, PDBDEV_00000020_resize, PDBDEV_00000021_resize, PDBDEV_00000022_resize, PDBDEV_00000023_resize, PDBDEV_00000024_resize, PDBDEV_00000025_resize, PDBDEV_00000026_resize, PDBDEV_00000027_resize, PDBDEV_00000028_resize, PDBDEV_00000029_resize, PDBDEV_00000031_resize, PDBDEV_00000032_resize, PDBDEV_00000033_resize"
Field_Separator=$IFS

IFS=", "
newtag="_1.png"
for filename in $mylist; do
    file="$filename.png"
    newfile="$filename$newtag"
    echo $file $newfile
    convert $file  +profile "*" -fuzz 1% -trim +repage -resize 500X500 -background white -gravity center -extent 500x500 $newfile

done
