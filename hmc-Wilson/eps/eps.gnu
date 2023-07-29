# adjust_eps() updates eps every 30 measurement steps
# Hip / 2023-07-29

set output
set terminal pdfcairo #size 4in, 3in
set output "eps.pdf"

set title "HMC-Wilson / 2-flavor / acceptance rate 80%"

set xlabel "measurement"
set ylabel "eps"

plot [0:1000][0.0:0.022] \
    "b02000k25000-eps.txt" w l, \
    "b04000k25000-eps.txt" w l, \
    "b06000k25000-eps.txt" w l

quit