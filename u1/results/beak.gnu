set output
set terminal pdfcairo #size 4in, 3in
set output "fig-beak-diagram.pdf"

set xlabel "1 / {/Symbol b}"
set ylabel "{/Symbol b} {/Symbol c}_T"

set key right center

plot "b1chiTsin.lst" w l dt "." lc "blue" title "analytic (Q_S)", \
  "b1chiTsin.col" w yerr pt 4 lc "blue" title "measured (Q_S)", \
  "b1chiT.lst" w l dt "." lc "red" title "numerical (Q_T)", \
  "b1chiT.col" w yerr pt 6 lc "red" title "measured (Q_T)", \
   "<echo 0 0.02533" w p pt 2 lc "magenta" title "continuum"

quit
