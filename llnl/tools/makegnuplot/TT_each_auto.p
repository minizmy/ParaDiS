# Gnuplot script file for plotting data in file "Torsion_output"
set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set y2tic auto                         # set ytics automatically

set title "Normalized Torque-Surface strain under torsion" font "Times-Roman,30"
unset colorbox

#set terminal pngcairo enhanced color font 'Times-Roman,10'
#set output 'Torque_den_Torsion.png'

set xlabel "Surface Strain" font "Times-Roman,30"  
set ylabel "Normalized Stress" font "Times-Roman,30"  
set y2label "Dislocation density" font "Times-Roman,30"   

set autoscale y
set autoscale y2

plot "Torsion_output" u 7:8 w lp axes x1y1 title "stress" , \
     "Torsion_output" u 7:9 w lp axes x1y2 title "density" 
pause 2
reread
