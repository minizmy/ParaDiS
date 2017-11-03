# Gnuplot script file for plotting data in file "stress_Total_strain"

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set y2tic auto                         # set ytics automatically
set title "Stress-strain under tension" 
unset colorbox

#set terminal png
#set output 'SS_Den_tension.png'

set xlabel "Strain"
set ylabel "Stress"
set y2label "Dislocation density"

set autoscale y
set autoscale y2

plot "stress_Total_strain" u 1:2 w lp axes x1y1 title "stress" , \
     "density" u 2:3 w lp axes x1y2 title "density" 
pause 5
reread
