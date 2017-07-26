# Gnuplot script file for plotting data in file "BEND_total"
set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Normalized Moment under pure bending" font "Times-Roman,15"
unset colorbox

set terminal pngcairo enhanced color font 'Times-Roman,10'
set output 'BEND.png'

set xlabel "K*R" font "Times-Roman,15"  
set ylabel "Normalized Moment (M/D^3)" font "Times-Roman,15"  

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
plot "D150_1" u 7:8 ls 1 title "D=150nm, 1", \
     "D150_2" u 7:8 ls 2 title "D=150nm, 2", \
     "D150_3" u 7:8 ls 2 title "D=150nm, 3", \
     "D150_4" u 7:8 ls 4 title "D=150nm, 4", \
     "D150_5" u 7:8 ls 5 title "D=150nm, 5"
pause 1
