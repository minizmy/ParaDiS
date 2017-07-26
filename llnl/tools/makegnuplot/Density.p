# Gnuplot script file for plotting data in file "TT_total"
set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "BCC under torsion" font "Times-Roman,15"
unset colorbox

set terminal pngcairo enhanced color font 'Times-Roman,10'
set output 'FT3_UN_total.png'

set xlabel "Surface Strain" font "Times-Roman,15"  
set ylabel "Normalized Stress" font "Times-Roman,15"  

#set xr[0:0.2]
#set yr[0:2e9]
#
# define line styles using explicit rgbcolor names
#
set style line 1 lt 2 lc rgb "red" lw 2 ps 0.2  pt 2
set style line 2 lt 2 lc rgb "blue" lw 2 ps 0.2 pt 2
set style line 3 lt 2 lc rgb "cyan" lw 2 ps 0.2 pt 2
set style line 4 lt 2 lc rgb "green" lw 2 ps 0.2 pt 2
set style line 5 lt 2 lc rgb "black" lw 2 ps 0.2 pt 2
#
show style line
#
plot "D1000_LT7_UN3_1" u 6:7 ls 1 title "D = 1000nm", \
     "D1000_LT7_UN3_2" u 6:7 ls 1 not, \
     "D1000_LT7_UN3_3" u 6:7 ls 1 not, \
     "D1000_LT7_UN3_4" u 6:7 ls 1 not, \
     "D1000_LT7_UN3_5" u 6:7 ls 1 not, \
     "D300_LT7_UN3_1" u 6:7 ls 2 title "D = 300nm", \
     "D300_LT7_UN3_2" u 6:7 ls 2 not, \
     "D300_LT7_UN3_3" u 6:7 ls 2 not, \
     "D300_LT7_UN3_4" u 6:7 ls 2 not, \
     "D300_LT7_UN3_5" u 6:7 ls 2 not, \
     "D600_LT7_UN3_1" u 6:7 ls 4 title "D = 600nm", \
     "D600_LT7_UN3_2" u 6:7 ls 4 not, \
     "D600_LT7_UN3_3" u 6:7 ls 4 not, \
     "D600_LT7_UN3_4" u 6:7 ls 4 not, \
     "D600_LT7_UN3_5" u 6:7 ls 4 not, \
     "D1000_LT7_UN3_1" u 6:7 ls 5 title "D = 1000nm", \
     "D1000_LT7_UN3_2" u 6:7 ls 5 not, \
     "D1000_LT7_UN3_3" u 6:7 ls 5 not, \
     "D1000_LT7_UN3_4" u 6:7 ls 5 not, \
     "D1000_LT7_UN3_5" u 6:7 ls 5 not
pause 1
