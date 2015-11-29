#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Provide the name of the text file with true values!" 
  exit 1
fi

./getResults.sh $1
make
./Main_plots.exe $1
python makeTexFiles.py $1
latex report.tex
dvips report.dvi 
ps2pdf report.ps
latex report_roofit.tex
dvips report_roofit.dvi 
ps2pdf report_roofit.ps
