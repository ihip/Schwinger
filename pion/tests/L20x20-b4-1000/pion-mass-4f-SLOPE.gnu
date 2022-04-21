set terminal pdfcairo
set size ratio 0.75
set output 'pion-L20x20-b4-1000-j1-chi2-pl4-jk20-4f-SLOPE.pdf'

set title "pion mass / j1 / L20x20 / {/Symbol b}=4 / 1000 meas."

set xlabel 'm'
set ylabel 'm_{/Symbol p}'

set key right bottom

# input parameters for small m Hosotani prediction
beta = 4.0
L = 20 # N_t

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

set arrow from b,0.0 to b,0.071 nohead lc rgb 'black' dt 2

plot [0:0.012][0:0.071] \
  k * x title "k_{2f}" lc "green", \
  k * x / 2.0 lc "cyan" title "k_{4f} = k_{2f} / 2", \
  k * x / sqrt(3.0) lc "cyan" dt 2 title "k_{4f} = k_{2f} / sqrt(3)", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-2f-HOR.col" u 1:2:3 w yerr pt 13 lc "red" ps 0.4 lw 0.3 \
    title "reweighted hyp. overlap (N_f=2)", \
  "pion-L20x20-b4-1000-j1-chi2-pl4-jk20-4f-HOR.col" u 1:2:3 w yerr pt 7 lc "black" ps 0.3 lw 0.3 \
    title "reweighted hyp. overlap (N_f=4)"

quit
