# Gnuplot script file for plotting data in file "filename.png"
# usage : >> gnuplot -e "filename='Torsion_output'" TT_each_filename.p
#
set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set y2tic auto                         # set ytics automatically

set title "Normalized Torque-Surface total strain" font "Times-Roman,15"
#set title filename
#set title sprintf('%s',filename)

set terminal pngcairo enhanced color font 'Times-Roman,10'
set output sprintf('%s.png',filename)

set xlabel "Surface Strain" font "Times-Roman,15"  
set ylabel "Normalized Stress" font "Times-Roman,15"  
set y2label "Dislocation density" font "Times-Roman,15"   

set autoscale y
set autoscale y2

plot filename u 7:8 w lp axes x1y1 title "stress" , \
     filename u 7:9 w lp axes x1y2 title "density" 
pause 5
