#!/bin/bash

# move files around so paraview will recognize them as a timeseries

for i in `seq 1 100000`
do
  if [[ -d mesh$i ]]
  then
    echo "directory exists"
    cp -v "./mesh$i/mesh$i.pvtu" "./mesh$i.pvtu"
    sed -i "s/0\/0.vtu/mesh$i\/0\/0.vtu/g" "./mesh$i.pvtu"
  fi
done

 #<Piece Source="0/0.vtu"/>

