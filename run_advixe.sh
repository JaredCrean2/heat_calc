#!/bin/bash

output_dir="./advixe_results"
exe="./src/simple_house2"

rm -r $output_dir

advixe-cl -collect survey -project-dir $output_dir -- $exe > fout_advise
#advixe-cl -collect tripcounts -flop -project-dir $output_dir -- $exe > fout_advise
#advixe-cl -collect map -flop -project-dir $output_dir -- $exe > fout_advise
advixe-cl -collect roofline -project-dir $output_dir -- $exe > fout_advise

