# this script generates a plot for the degree wrt z axis for different data files


set encoding iso_8859_1 
#set title 'H-bonding Network Cumulative Degree' font "Times-Roman, 36"
#set xlabel "z-axis (\305)" font "Times-Roman, 28"
#set ylabel 'Degree' font "Times-Roman, 28"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xrange[0:14.41]
set yrange[2:5]
plot 'PBE/degree_z.dat' title 'PBE' w lines lw 3 linecolor rgb "red", 'LDA/degree_z.dat' title 'LDA' w lines lw 3 linecolor rgb "royalblue", 'vdW/degree_z.dat' title 'vdW-DF' w lines lw 3 linecolor rgb "olive", 'DFTB/degree_z.dat' title 'DFTB' w lines lw 3 linecolor rgb "brown"