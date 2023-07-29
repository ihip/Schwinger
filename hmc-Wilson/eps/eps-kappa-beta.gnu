# integration step (eps) dependence on kappa and beta
# Hip / 2023-07-29

set output
set terminal pdfcairo #size 4in, 3in
set output "eps-kappa-beta.pdf"

set title "HMC-Wilson / 2 flavor / nmeas = 1000"

set xlabel "kappa"
set ylabel "eps"

plot [0.2:0.3] \
    "beta2-eps.txt" u 1:2:3 w yerr pt 12 ps 0.5 lw 0.5, \
    "beta4-eps.txt" u 1:2:3 w yerr pt 12 ps 0.5 lw 0.5, \
    "beta6-eps.txt" u 1:2:3 w yerr pt 12 ps 0.5 lw 0.5

quit