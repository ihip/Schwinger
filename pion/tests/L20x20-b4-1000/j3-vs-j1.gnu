set terminal pdfcairo
set size ratio 0.75
set output 'j3-vs-j1-L20x20-b4-1000-chi2-pl4-jk20.pdf'

set title "pion mass / j3 vs. j1 / L20x20 / {/Symbol b}=4 / 1000 meas."

set xlabel 'm'
set ylabel 'm_{/Symbol p}'

set key right bottom

# input parameters for small m Hosotani prediction
beta = 4.0
L = 20 # N_t

g = 1 / sqrt(beta)
gamma = 0.5772156649 # Euler-Mascheroni constant
Meta = sqrt(2.0 / (pi * beta)) # eta mass

# bellow this boundary small m Hosotani prediction should be valid
b = 1.0 / (2.0 * L * sqrt(Meta * L))
print "b = ", b
# L = 20 & beta = 4.0 -> b = 0.00885054425344672

# slope of Hosotani prediction for small m
k = 4 * sqrt(2.0) * sqrt(Meta * L * exp(gamma) / (4 * pi))
print "k = ", k
# L = 20 & beta = 4.0 -> k = 6.01562668203976

plot [0:0.16][0:0.51] k * x title "Hosotani (small m)" lc "gray", \
  2.16333333 * (x**2 * g)**(1.0/3.0) title "Hosotani (large m)" lc "violet", \
  2.008 * (x**2 * g)**(1.0/3.0) title "Sine-Gordon" lc "light-blue", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-2f-WD.col" w yerr pt 7 lc "blue" ps 0.4 lw 0.3 \
    title "dynamical Wilson (j1, N_f=2)", \
  "pion-L20x20-b4-1000-j3-chi2-pl4-jk20-2f-WD.col" w yerr pt 6 lc "blue" ps 0.4 lw 0.3 \
    title "dynamical Wilson (j3, N_f=2)", \
  "pion-L20x20-b4-1000-chi2-pl4-jk20-2f-HOR.col" i 0 u 1:2:3 w yerr pt 13 lc "red" ps 0.4 lw 0.3 \
    title "reweighted hyp. overlap (j1, N_f=2)", \
  "pion-L20x20-b4-1000-chi2-pl4-jk20-2f-HOR.col" i 1 u 1:2:3 w yerr pt 12 lc "red" ps 0.4 lw 0.3 \
    title "reweighted hyp. overlap (j3, N_f=2)"

quit
