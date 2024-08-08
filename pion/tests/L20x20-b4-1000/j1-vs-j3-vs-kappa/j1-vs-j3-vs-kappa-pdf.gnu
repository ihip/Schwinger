set terminal pdfcairo
set size ratio 0.6
set output 'j1-vs-j3-vs-kappa-L20x20-b4-1000.pdf'

set title "pion mass / j1 vs. j3 / L20x20 / {/Symbol b}=4 / 1000 meas."

set xlabel '{/Symbol k}'
set ylabel 'm_{/Symbol p}'

set key left bottom

plot [0.198:0.267][0.0:1.2] \
  "mass-L20x20-b4-WD.col" i 0 u 1:2:3 w yerr pt 6 ps 0.4 lw 0.5 title columnheader, \
  "mass-L20x20-b4-WD.col" i 1 u 1:2:3 w yerr pt 12 ps 0.4 lw 0.5 title columnheader, \
  "mass-L20x20-b4-WD.col" i 2 u 1:2:3 w yerr pt 3 ps 0.4 lw 0.5 title columnheader(1)  

quit
