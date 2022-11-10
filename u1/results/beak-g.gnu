set output
set terminal pdfcairo #size 4in, 3in
# set term pdfcairo font "Times-New-Roman,12"
set output "BeakDiagram-g.pdf"

set xlabel "g"
set ylabel "{/Symbol c}_t^q / g^2"

set key right center

plot "b1chiTsin.lst" u (sqrt($1)):2 w l dt "." lc "blue" title "analytic (Q_S)", \
  "b1chiTsin.col" u (sqrt($1)):2 w yerr pt 4 lc "blue" title "simulation (Q_S)", \
  "b1chiT.lst" u (sqrt($1)):2 w l dt "." lc "red" title "analytic (Q_T)", \
  "b1chiT.col" u (sqrt($1)):2 w yerr pt 6 lc "red" title "simulation (Q_T)", \
  "<echo 0 0.02533" w p pt 2 lc "magenta" title "continuum"

quit
