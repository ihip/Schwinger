set terminal pdfcairo
set size ratio 0.75
set output 'sigma-L20x20-b4-1000-jk25.pdf'

set xlabel 'm / g'
set ylabel '{/Symbol S} / g'

set key right bottom

plot [0:0.3][0:0.4] \
  "b4-Nf0-1000-jk25.sig" w yerr pt 7 lc "red" ps 0.3 lw 0.3 title "N_f = 0", \
  "b4-Nf1-1000-jk25.sig" w yerr pt 7 lc "black" ps 0.3 lw 0.3 title "N_f = 1", \
  "b4-Nf2-1000-jk25.sig" w yerr pt 7 lc "blue" ps 0.3 lw 0.3 title "N_f = 2"

quit
