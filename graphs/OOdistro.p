# this script generates a plot for the nearest neighbor OO bond length for different data files


set encoding iso_8859_1 
#set title 'O-O Nearest Neighbour Distance Distribution' font "Times-Roman, 36"
#set xlabel "Distance (\305)" font "Times-Roman, 24"
#set ylabel "Probability (1/\305)" font "Times-Roman, 24" offset -2
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xrange[2.2:3.5]
plot 'PBE/OOdistro.dat' title 'PBE' w lines lw 3 linecolor rgb "red", 'LDA/OOdistro.dat' title 'LDA' w lines lw 3 linecolor rgb "royalblue", 'vdW/OOdistro.dat' title 'vdW-DF' w lines lw 3 linecolor rgb "olive", 'DFTB/OOdistro.dat' title 'DFTB' w lines lw 3 linecolor rgb "brown"