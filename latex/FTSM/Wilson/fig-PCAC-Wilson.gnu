
set terminal pdfcairo size 4in, 3in

# set aspect ratio to 1
#set size ratio -1

# set output file name
set output "fig-PCAC-Wilson.pdf"

# set title "PCAC"
set xlabel "{/Symbol k}"
set ylabel "m"

plot [0.2:0.3] \
"fig-PCAC-Wilson.dat" i 0 u 1:2:3 w yerr lc "blue" pt 6 ps 0.2 title columnhead, \
"fig-PCAC-Wilson.dat" i 1 u 1:2:3 w yerr lc "violet" pt 12 ps 0.2 title columnhead, \
0 lc "gray" title ""

quit

