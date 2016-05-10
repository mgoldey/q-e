#!/usr/bin/gnuplot
reset
# wxt
set terminal wxt size 1000,750 enhanced font 'Times-Roman,28' persist
# png
#set terminal pngcairo size 1000,750 enhanced font 'Times-New-Roman,28'
#set output 'bc3_zn_coupling.png'
#set terminal postscript eps size 6.0, 4.50 enhanced  color font 'Times-New-Roman,28'
#set output 'vtot.eps'
# svg
#set terminal epslatex size 800,600 font 'Times-Roman, 18'
#set output 'sexc.eps'

# define axis
# remove border on top and right 
set style line 11 lc rgb '#000000' lt 1
set border 3 back ls 11
set tics nomirror
# define grid
#set style line 12 lc rgb '#000000' lt 0 lw 1
#set grid back ls 12

# color definitions
set style line 1 lc rgb '#00008B' pt 3 ps 2.2 lt 1 lw 2
set style line 2 lc rgb 'red' pt 7 ps 2.2 lt 1 lw 2
set style line 3 lc rgb 'blue' pt 9 ps 2.2 lt 1 lw 2
set style line 4 lc rgb 'black' pt 11 ps 2.2 lt 1 lw 2
set style line 5 lc rgb '#006400' pt 13 ps 2.2 lt 1 lw 2
set style line 6 lc rgb '#708090' pt 5 ps 2.2 lt 1 lw 2

set pointsize 1.5


#set key invert at 9.5,16.2
set key bottom right

set xlabel 'd [Bohr]'
set title 'Surface Interaction Energy'
set ylabel 'Energy [Ry]'
#unset ytics
#set xrange [2.9:17]
#set yrange [2.9:17]
#set label 'b)' at 6.80, 3.2 
#set label 'a)' at 6.75, 2.1 

#set label "Mean \% error \n------------------\n SEXC_o 2.6\% \n G_oW_o at SEXC_o  3.5\% \n G_oW_o at PBE  4.0\%" enhanced font 'Times-Roman,16' at 6,16
#set label "   \% Error \n------------------" at 13,8

#f1(x) = a1*exp(-b1*x/2.0) + c1 
# a1 = 1.0; b1 = 1.0; c1 = 0.1 
# f2(x) = a2*exp(-b2*x/2.0) + c2 
# a2 = 1.0; b2 = 1.0; c2 = 0.1 

#fit f1(x) 'bc3' using 1:2 via a1, b1, c1
# fit f2(x) 'bc1' using 1:2 via a2, b2, c2

e2 = 2.0 #Ry units

#plot 'bc3' u 1:2 w p ls 1 t 'bc3', a1*exp(-b1*x/2.0)+c1 notitle, 'bc1' u 1:2 w p ls 2 t 'bc1', a2*exp(-b2*x/2.0)+c2 notitle,
plot 'data' u 1:2 w p ls 1 notitle , -e2/(2*x) ls 1 t '-e^2/(2d)'
