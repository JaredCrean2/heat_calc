#!/bin/bash

dest=$1

mkdir -vp $dest
mv -v ./mesh* $dest/
cp -v ./fout $dest/
cp -v ./simple_house_data.txt $dest/
cp -v ./*.png $dest/
