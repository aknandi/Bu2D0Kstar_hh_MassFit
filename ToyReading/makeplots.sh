./getResults.sh $1 $2
make
./Main_plots.exe $1 $2
latex report.tex
dvips report.dvi -o report.ps
ps2pdf report.ps
