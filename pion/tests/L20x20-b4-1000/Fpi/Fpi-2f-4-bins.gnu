set terminal pdfcairo
set size ratio 0.75
set output 'Fpi-2f-4-bins.pdf'

set title "F_{/Symbol p} / 2 flavors / j1 / L20x20 / {/Symbol b}=4 / 4 x 250 meas."

set xlabel 'm'
set ylabel 'F_{/Symbol p}'

set key right bottom

plot [0:0.51][0:0.61] \
  "bin1-1-250.col" i 0 u 1:6:7 w yerr pt 13 lc "cyan" ps 0.4 lw 0.3 \
     title "bin 1", \
  "bin2-251-500.col"i 0 u 1:6:7 w yerr pt 13 lc "red" ps 0.4 lw 0.3 \
     title "bin 2", \
  "bin3-501-750.col" i 0 u 1:6:7 w yerr pt 13 lc "green" ps 0.4 lw 0.3 \
     title "bin 3", \
  "bin4-751-1000.col" i 0 u 1:6:7 w yerr pt 13 lc "orange" ps 0.4 lw 0.3 \
     title "bin 4", \
  "../pion-L20x20-b4-1000-j1-chi2-pl4-jk20-2f-HOR.col" u 1:6:7 w yerr \
     pt 7 lc "black" ps 0.4 lw 1.0 title "all together"

quit
