# Gnuplot script file for plotting data in file "Torsion_output"

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Torsion, D=1000nm"
unset colorbox

set terminal png
set output 'TP_D1000.png'

set xlabel "Surface Plastic Strain"
set ylabel "Normalized Stress"

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
plot "Torsion_D1000_1" u 6:8 ls 1 title "case 1", \
     "Torsion_D1000_2" u 6:8 ls 2 title "case 2", \
     "Torsion_D1000_3" u 6:8 ls 3 title "case 3", \
     "Torsion_D1000_4" u 6:8 ls 4 title "case 4", \
     "Torsion_D1000_5" u 6:8 ls 5 title "case 5"
pause 1
