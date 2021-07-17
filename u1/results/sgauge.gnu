set output
set terminal win

plot "sgauge.out" u 1:2:3 w yerr pt 6

pause -1
quit
