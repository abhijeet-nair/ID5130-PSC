set terminal postscript eps enhanced "Times" 25

set key bottom right 
set key spacing 1
set key samplen 1

set ylabel 'Time (s)'
set xlabel 'Matrix size N'

set size 1.0,1.0
set xtics
set ytics 

set output 'Matrix_Product_Time.eps

plot 'prod.txt' using 1:2 with linespoints lt 1 pt 5 ps 3 title 'Matrix product'
