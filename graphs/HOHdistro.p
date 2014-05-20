# this script generates a plot for the H-O-H angle for different data files


set encoding iso_8859_1 
#set title 'H-O-H Bending Angle Distribution' font "Times-Roman, 36"
#set xlabel "Angle (\260)" font "Times-Roman, 32"
#set ylabel "Probability (1/\260)" font "Times-Roman, 32"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xrange[80:140]
plot 'PBE/HOHdistro.dat' title 'PBE' w lines lw 3 linecolor rgb "red", 'LDA/HOHdistro.dat' title 'LDA' w lines lw 3 linecolor rgb "royalblue", 'vdW/HOHdistro.dat' title 'vdW-DF' w lines lw 3 linecolor rgb "olive", 'DFTB/HOHdistro.dat' title 'DFTB' w lines lw 3 linecolor rgb "brown"