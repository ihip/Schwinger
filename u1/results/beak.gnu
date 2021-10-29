set output
set terminal pdfcairo #size 4in, 3in
set output "BeakDiagram.pdf"

set xlabel "1 / {/Symbol b}"
set ylabel "{/Symbol b} {/Symbol c}_t^q"

set key right center

plot "b1chiTsin.lst" w l dt "." lc "blue" title "analytic (Q_S)", \
  "b1chiTsin.col" w yerr pt 4 lc "blue" title "simulation (Q_S)", \
  "b1chiT.lst" w l dt "." lc "red" title "analytic (Q_T)", \
  "b1chiT.col" w yerr pt 6 lc "red" title "simulation (Q_T)", \
  "<echo 0 0.02533" w p pt 2 lc "magenta" title "continuum"

quit
