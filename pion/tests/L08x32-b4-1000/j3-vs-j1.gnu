set terminal pdfcairo
set size ratio 0.75
set output 'j3-vs-j1-L08x32-b4-1000-chi2-pl4-jk20.pdf'

set title "pion mass / j3 vs. j1 / L08x32 / {/Symbol b}=4 / 1000 meas."

set xlabel 'm'
set ylabel 'm_{/Symbol p}'

set key right bottom

plot [0:0.31][0:0.91] 2.16333333 * (x**2 / 2.0)**(1.0/3.0) title "Hosotani", \
  "pion-L08x32-b4-1000-j1-chi2-pl4-jk20-2f-WD.col" w yerr pt 7 lc "blue" ps 0.4 lw 0.3 \
    title "dynamical Wilson (j1, N_f=2)", \
  "pion-L08x32-b4-1000-j3-chi2-pl4-jk20-2f-WD.col" w yerr pt 6 lc "blue" ps 0.4 lw 0.3 \
    title "dynamical Wilson (j3, N_f=2)", \
  "pion-L08x32-b4-1000-chi2-pl4-jk20-2f-HOR.col" i 0 u 1:2:3 w yerr pt 13 lc "red" ps 0.4 lw 0.3 \
    title "reweighted hyp. overlap (j1, N_f=2)", \
  "pion-L08x32-b4-1000-chi2-pl4-jk20-2f-HOR.col" i 1 u 1:2:3 w yerr pt 12 lc "red" ps 0.4 lw 0.3 \
    title "reweighted hyp. overlap (j3, N_f=2)"

quit
