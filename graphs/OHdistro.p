# this script generates a plot for OH bond length for different data files


set encoding iso_8859_1 
#set title 'O-H Bond Length Distribution' font "Times-Roman, 36"
#set xlabel "Distance (\305)" font "Times-Roman, 32"
#set ylabel "Probability (1/\305)" font "Times-Roman, 32"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xrange[0.89:1.21]
plot 'PBE/OHdistro.dat' title 'PBE' w lines lw 3 linecolor rgb "red", 'LDA/OHdistro.dat' title 'LDA' w lines lw 3 linecolor rgb "royalblue", 'vdW/OHdistro.dat' title 'vdW-DF' w lines lw 3 linecolor rgb "olive", 'DFTB/OHdistro.dat' title 'DFTB' w lines lw 3 linecolor rgb "brown"