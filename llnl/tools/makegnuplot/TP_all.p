# Gnuplot script file for plotting data in file "Torsion_output"
#set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "FCC under torsion"
unset colorbox

set terminal png
set output 'TP_all.png'

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
plot "Torsion_D150_1" u 6:8 ls 1 title "D = 150nm", \
     "Torsion_D150_2" u 6:8 ls 1 not, \
     "Torsion_D150_3" u 6:8 ls 1 not, \
     "Torsion_D150_4" u 6:8 ls 1 not, \
     "Torsion_D150_5" u 6:8 ls 1 not, \
     "Torsion_D300_1" u 6:8 ls 2 title "D = 300nm", \
     "Torsion_D300_2" u 6:8 ls 2 not, \
     "Torsion_D300_3" u 6:8 ls 2 not, \
     "Torsion_D300_4" u 6:8 ls 2 not, \
     "Torsion_D300_5" u 6:8 ls 2 not, \
     "Torsion_D600_1" u 6:8 ls 4 title "D = 600nm", \
     "Torsion_D600_2" u 6:8 ls 4 not, \
     "Torsion_D600_3" u 6:8 ls 4 not, \
     "Torsion_D600_4" u 6:8 ls 4 not, \
     "Torsion_D600_5" u 6:8 ls 4 not, \
     "Torsion_D1000_2" u 6:8 ls 5 title "D = 1000nm", \
     "Torsion_D1000_3" u 6:8 ls 5 not, \
     "Torsion_D1000_4" u 6:8 ls 5 not, \
     "Torsion_D1000_5" u 6:8 ls 5 not
pause 1
