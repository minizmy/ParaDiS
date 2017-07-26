# Gnuplot script file for plotting data in file "Ts_Ps_SS_Den"

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Tension, D=150nm"
unset colorbox

set terminal png
set output 'SS_D150.png'

#set xlabel "Strain" font "Times-Roman,15"  
#set ylabel "Stress" font "Times-Roman,15"  

set xlabel "Strain"
set ylabel "Stress"

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
plot "Tension_D150_1" u 1:3 ls 1 title "case 1", \
     "Tension_D150_2" u 1:3 ls 2 title "case 2", \
     "Tension_D150_3" u 1:3 ls 3 title "case 3", \
     "Tension_D150_4" u 1:3 ls 4 title "case 4", \
     "Tension_D150_5" u 1:3 ls 5 title "case 5"
pause 1
