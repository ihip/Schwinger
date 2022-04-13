set terminal pdfcairo
set size ratio 0.75
set output 'pion-L20x20-b4-1000-j1-chi2-pl4-jk20.pdf'

set title "pion mass / j1 / L20x20 / {/Symbol b}=4 / 1000 meas."

set xlabel 'm'
set ylabel 'm_{/Symbol p}'

set key right bottom

plot [0:0.51][0:1.3] 2.16333333 * (x**2 / 2.0)**(1.0/3.0) title "Hosotani", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-QM.dat" w yerr pt 7 lc "orange" ps 0.3 lw 0.3 title "quenched", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-WD.col" w yerr pt 7 lc "blue" ps 0.2 lw 0.3 title "dynamical Wilson (N_f=2)", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-HOR.col" w yerr pt 13 lc "red" ps 0.4 lw 0.3 title "reweighted hypercube overlap (N_f=2)"

quit
