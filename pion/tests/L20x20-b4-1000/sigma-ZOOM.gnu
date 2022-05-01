set terminal pdfcairo
set size ratio 0.75
set output 'sigma-L20x20-b4-1000-jk20-ZOOM.pdf'

set title "chiral condensate / L20x20 / {/Symbol b}=4 / 1000 meas."

set xlabel 'm'
set ylabel '{/Symbol S}'

set key right bottom

# input parameters for small m Hosotani prediction
beta = 4.0
L = 20 # N_t

gamma = 0.5772156649 # Euler-Mascheroni constant
mu = sqrt(2.0 / (pi * beta)) # eta mass

# bellow this boundary small m Hosotani prediction for Nf = 2 should be valid
b = 1.0 / (L * sqrt(mu * L))
print "b = ", b
#L = 20 & beta = 4.0 -> b = 0.0177010885068934

# slope of Hosotani prediction for small m (2 flavors, PLB 350, eq. 36)
k2 = 2 * exp(gamma) * mu * L / pi**2

# slope of Hosotani prediction for small m (Nf flavors, PRD 53, eq. 16)
k(Nf) = 2 * Nf * (mu * L * exp(gamma) / (4 * pi))**(2.0 / Nf) / (pi * (Nf - 1))
print "PLB 350 N_f = 2: k = ", k2
print "PRD 53: k(N_f = 3) = ", k(3)
print "PRD 53: k(N_f = 4) = ", k(4)

set arrow from b,0.0 to b,0.071 nohead lc rgb 'black' dt 2

plot [0:0.011][0:0.011] \
  k2 * x lc "red" dt 2 title "PLB 350: N_f = 2", \
  k(3.0) * x lc "blue" dt 2 title "PRD 53: N_f = 3", \
  k(4.0) * x lc "green" dt 2 title "PRD 53: N_f = 4", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-2f-HOR.col" u 1:4:5 w yerr pt 13 lc "red" ps 0.4 lw 0.3 \
     title "reweighted HO (N_f = 2)", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-3f-HOR.col" u 1:4:5 w yerr pt 7 lc "blue" ps 0.3 lw 0.3 \
     title "reweighted HO (N_f = 3)", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-4f-HOR.col" u 1:4:5 w yerr pt 13 lc "green" ps 0.4 lw 0.3 \
     title "reweighted HO (N_f = 4)"

quit
