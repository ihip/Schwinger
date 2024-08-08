set output
set terminal win

plot [0.2:0.27][0.0:1.2] \
  "mass-L20x20-b4-WD.col" i 0 u 1:2:3 w yerr pt 6 title columnheader, \
  "mass-L20x20-b4-WD.col" i 1 u 1:2:3 w yerr pt 12 title columnheader(1), \
  "mass-L20x20-b4-WD.col" i 2 u 1:2:3 w yerr pt 3 title columnheader(1)  

pause -1
quit