set output
set terminal pdfcairo #size 4in, 3in
set output "acceptance.pdf"

set xlabel "iteration"
set ylabel "acceptance rate"

set key right bottom
# unset key

plot "b02k025eps1.acc" w l, \
    "b04k025eps1.acc" w l, \
    "b06k025eps1.acc" w l    

quit