# this script generates a plot for the O-O pair correlation function for different data files


set encoding iso_8859_1 
#set title 'O-O Pair' font "Times-Roman, 36"
#set xlabel "Distance (\305)" font "Times-Roman, 24"
#set ylabel "Probability (1/\305)" font "Times-Roman, 24" offset -2
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xrange[0:5.0]
set yrange[-5:*]
plot 'PBE/OH_RDF.dat' title 'PBE' w lines lw 3 linecolor rgb "red", 'LDA/OH_RDF.dat' title 'LDA' w lines lw 3 linecolor rgb "royalblue", 'vdW/OH_RDF.dat' title 'vdW-DF' w lines lw 3 linecolor rgb "olive", 'DFTB/OH_RDF.dat' title 'DFTB' w lines lw 3 linecolor rgb "brown"
