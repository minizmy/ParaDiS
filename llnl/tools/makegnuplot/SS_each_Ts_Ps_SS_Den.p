# Gnuplot script file for plotting data in file "Ts_Ps_TS_Den"
# usage : >> gnuplot SS_each_Ts_Ps_SS_Den.p
#

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set y2tic auto                         # set ytics automatically

set title "Stress-Strain under Tension"
#set title filename
#set title sprintf('%s',filename)

set terminal png
set output sprintf('%s.png',"Ts_Ps_SS_Den2")

set xlabel "Strain"
set ylabel "Stress"
set y2label "Dislocation density"

set autoscale y
set autoscale y2

plot "Ts_Ps_SS_Den" u 1:3 w lp axes x1y1 title "stress" , \
     "Ts_Ps_SS_Den" u 1:4 w lp axes x1y2 title "density"
pause 5
