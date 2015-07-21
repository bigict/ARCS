#!/bin/bash
outf=$1
path=$2
title=$3
echo " set terminal postscript color eps 
 #set multiplot
 set xlabel 'reference pos(bp)'
 set ylabel 'Number of coverage'
 set xrange [-101:700]
 set yrange [0:50]
 set term pdfcairo lw 2 font 'Times New Roman,8'
 set size 0.9,0.8
 set grid
 set output '$1'
 plot '$2' u 1:3 with linespoints linetype 1 linecolor 1 linewidth 0.5 pointtype 9 pointsize 0.5 smooth csplines  title '010', '$2'  u 1:5 with linespoints linetype 1 linecolor 5 linewidth 0.5 pointtype 9 pointsize 0.5 smooth csplines  title '030'	,'$2' u 1:7 with linespoints linetype 1 linecolor 0 linewidth 0.5 pointtype 9 pointsize 0.5 smooth csplines  title '000','$2' u 1:9 with linespoints linetype 1 linecolor 9 linewidth 0.5 pointtype 9 pointsize 0.5 smooth csplines  title '020'
" > temp
gnuplot temp
#, './33mer/repeat_count' u 1:2 with linespoints linetype 1 linecolor 1 linewidth 0.5 pointtype 7 pointsize 0.5 t 'K=33'
#set output './ecoli_points.eps'
#lot './23mer/repeat_count' u 1:2 with points linetype 1 linecolor 3 linewidth 0.5 pointtype 9 pointsize 0.5 t 'K=23', './33mer/repeat_count' u 1:2 with points linetype 1 linecolor 1 linewidth 0.5 pointtype 7 pointsize 0.5 t 'K=33'
# plot '$2' u 1:2 with linespoints linetype 1 linecolor 2 linewidth 0.5 pointtype 9 pointsize 0.5 smooth csplines  title '000' , u 1:3 with linespoints linetype 1 linecolor 3 linewidth 0.5 pointtype 9 pointsize 0.5 smooth csplines  title '010',  u 1:4 with linespoints linetype 1 linecolor 3 linewidth 0.5 pointtype 9 pointsize 0.5 smooth csplines  title '020'
