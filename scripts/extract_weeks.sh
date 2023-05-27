#!/bin/bash

months=("jan" "feb" "mar" "apr" "may" "jun" "july" "aug" "sep" "oct" "nov" "dec")

for i in `seq 0 11`
do
  mon="${months[$i]}"
  echo "month = $mon"
  i2=$(( i + 1 ))
  data_extractor ./simple_house_data.txt weather.wea "$i2/1/2010" "$i2/7/2010" simple_house_"$mon"_week1.txt
  data_extractor ./simple_house_data.txt weather.wea "$i2/8/2010" "$i2/14/2010" simple_house_"$mon"_week2.txt
  data_extractor ./simple_house_data.txt weather.wea "$i2/15/2010" "$i2/21/2010" simple_house_"$mon"_week3.txt
  data_extractor ./simple_house_data.txt weather.wea "$i2/22/2010" "$i2/28/2010" simple_house_"$mon"_week4.txt
done
