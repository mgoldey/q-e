#!/usr/bin/gnuplot
reset
# wxt
set terminal wxt size 800,600 enhanced font 'Times-Roman,16' persist
# png
#set terminal pngcairo size 1000,750 enhanced font 'Times-New-Roman,24'
#set output 'sexc.png'
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
set style line 1 lc rgb 'black' pt 6 ps 1 lt 1 lw 2 

set pointsize 1.5


set key  bottom left

set xlabel 'd [Ang]'
set ylabel '|Hab| [mHa]'

set logscale y 10
#unset ytics
set xrange [4.5:9.5]
set yrange [0.01:10]
#set label 'b)' at 6.80, 3.2 
#set label 'a)' at 6.75, 2.1 

#set label "Mean \% error \n------------------\n SEXC_o 2.6\% \n G_oW_o at SEXC_o  3.5\% \n G_oW_o at PBE  4.0\%" enhanced font 'Times-Roman,16' at 6,16
#set label "   \% Error \n------------------" at 13,8

# SEXC IS SEXCo at PBE
plot 'results' u 1:($2*1000) w p pt 7 lc rgb 'black' t 'CDFT QE' ,\
'WuVan' u 1:2 w p pt 7 lc rgb 'blue' t 'Wu, Van-Voorhis',\
'cpmd' u 1:2 w p pt 7 lc rgb 'red' t 'CPMD'
