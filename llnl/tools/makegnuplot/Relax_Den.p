# Gnuplot script file for plotting data in file "density"
set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Relaxation density" font "Times-Roman,15"
unset colorbox

#set terminal pngcairo enhanced color font 'Times-Roman,10'
#set output 'BCC_Tension_D1000.png'

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
plot "density" u 2:3 ls 1 
pause -1
