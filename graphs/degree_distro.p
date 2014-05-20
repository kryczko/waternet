# this script generates a plot for the degree distro for different data files


set encoding iso_8859_1 
#set title 'H-bonding Network Cumulative Degree Distribution' font "Times-Roman, 36"
#set xlabel "Degree" font "Times-Roman, 32"
#set ylabel 'Probability [1/Degree]' font "Times-Roman, 32"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xrange[0:8]
set yrange[0:1]
plot 'PBE/degree_distro.dat' every 2 title 'PBE' w lines lw 3 linecolor rgb "red", 'LDA/degree_distro.dat' every 2 title 'LDA' w lines lw 3 linecolor rgb "royalblue", 'vdW/degree_distro.dat' every 2 title 'vdW-DF' w lines lw 3 linecolor rgb "olive", 'DFTB/degree_distro.dat' every 2 title 'DFTB' w lines lw 3 linecolor rgb "brown"