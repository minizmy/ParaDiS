# Gnuplot script file for plotting data in file "Bending_output"
set termoption dash

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set y2tic auto                         # set ytics automatically

set title "Normalized Moment under pure bending" font "Times-Roman,15"
unset colorbox

#set terminal pngcairo enhanced color font 'Times-Roman,10'
#set output 'Bending.png'

set xlabel "K*R" font "Times-Roman,15"  
set ylabel "Normalized Moment (M/D^3)" font "Times-Roman,15"  
set y2label "Dislocation density" font "Times-Roman,15"   

set autoscale y
set autoscale y2

plot "Bending_output" u 7:8 w lp axes x1y1 title "stress" , \
     "Bending_output" u 7:9 w lp axes x1y2 title "density" 
pause 2
reread
