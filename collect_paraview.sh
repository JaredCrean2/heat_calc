#!/bin/bash

# make symbolic links so paraview recognizes the files
# as a sequence

for i in `seq 0 100`
do
  cp -vs ./mesh$i/mesh$i.pvtu ./mesh$i.pvtu
  sed -i "s/Piece Source=\"/Piece Source=\"mesh$i\//g" ./mesh$i.pvtu
done
