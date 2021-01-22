
set terminal pdfcairo size 4in, 3in

# set output file name
set output "Nt32b02k27-all.pdf"

set title "N_{t} = 32, {/Symbol b} = 2.0, {/Symbol k} = 0.27, nmeas = 1000"
set xlabel "L"
set ylabel "m_{{/Symbol p}}"

plot [:][0:2] \
"Nt32b02k27-all.dat" i 0 u 1:3:4 w yerr pt 6 ps 0.5 title columnhead, \
"Nt32b02k27-all.dat" i 1 u 1:3:4 w yerr pt 12 ps 0.5 title columnhead, \
"Nt32b02k27-all.dat" i 2 u 1:3:4 w yerr pt 2 ps 0.5 title columnhead, \
"Nt32b02k27-all.dat" i 3 u 1:3:4 w yerr pt 4 ps 0.5 title columnhead

quit

