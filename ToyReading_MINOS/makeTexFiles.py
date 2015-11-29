## make the necessary tex files to plot all variables
import sys

if len(sys.argv)<2:
  sys.exit("Provide an argument!")

## read in text file
infile = open(sys.argv[1])

# output file
sys.stdout = open('report.tex', 'w')
print '\documentclass[11pt,a4paper]{article}'
print '\usepackage{graphicx}'
print '\usepackage{makeidx}'
print '\usepackage{lscape}'
print '\pagestyle{plain}'
print '\pagenumbering{arabic}'
print '\\topmargin 0pt'
print '\headheight 12pt' 
print '\\footskip 40pt'
print '\\textheight 24cm'
print '\\textwidth 17.2cm'
print '\oddsidemargin 0mm'
print '\evensidemargin 0mm'
print '\parskip 3mm'
print '\parindent 0mm'
print '\makeindex'
print '\usepackage{amssymb}'
print '\usepackage{lscape}'
print '\\begin{document}'

prefix='\includegraphics[width=0.99\\textwidth]{results/'
suffix='_plots.eps}\n'
for line in infile:
  entry = line.split()
  outline = prefix+entry[0]+suffix
  print outline

print '\clearpage'
print '\includegraphics[width=0.5\\textwidth]{results/matrix.eps}'
print '\end{document}'

# output file for roofit plots
sys.stdout = open('report_roofit.tex', 'w')
print '\documentclass[11pt,a4paper]{article}'
print '\usepackage{graphicx}'
print '\usepackage{makeidx}'
print '\usepackage{lscape}'
print '\pagestyle{plain}'
print '\pagenumbering{arabic}'
print '\\topmargin 0pt'
print '\headheight 12pt' 
print '\\footskip 40pt'
print '\\textheight 24cm'
print '\\textwidth 17.2cm'
print '\oddsidemargin 0mm'
print '\evensidemargin 0mm'
print '\parskip 3mm'
print '\parindent 0mm'
print '\makeindex'
print '\usepackage{amssymb}'
print '\usepackage{lscape}'
print '\\begin{document}'

prefix='\includegraphics[width=0.45\\textwidth]{results_roofit/'
suffix='_plots.eps}'

infile.seek(0) # go back to start of file
counter=1
for line in infile:
  entry = line.split()
  outline = prefix+entry[0]+suffix
  if(counter%2==0):
    outline = outline+'\n'
  print outline
  counter=counter+1

print '\end{document}'






