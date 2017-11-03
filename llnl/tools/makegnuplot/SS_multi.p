# Gnuplot script file for plotting data in file "TT_total"
set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "BCC tension D=1000nm" font "Times-Roman,15"
unset colorbox

set terminal pngcairo enhanced color font 'Times-Roman,10'
set output 'BCC_Tension_D1000.png'

set xlabel "Strain" font "Times-Roman,15"  
set ylabel "Stress" font "Times-Roman,15"  

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
plot "D1000_LT7_1" u 1:2 ls 1 title "1", \
     "D1000_LT7_2" u 1:2 ls 2 title "2", \
     "D1000_LT7_3" u 1:2 ls 3 title "3", \
     "D1000_LT7_4" u 1:2 ls 4 title "4", \
     "D1000_LT7_5" u 1:2 ls 5 title "5"
pause 1
