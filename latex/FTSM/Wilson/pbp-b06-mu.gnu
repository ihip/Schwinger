
set terminal pdfcairo size 4in, 3in

# set aspect ratio to 1
#set size ratio -1

# set output file name
set output "pbp-b06-mu.pdf"

#set title "Wilson fermions, {/Symbol b} = 6.0"
set xlabel "T / {/Symbol m}"
set ylabel "<~{/Symbol Y}{.8-}{/Symbol Y}>"

set logscale x

# pt = pointtype (12 = diamonds); ps = pointsize
plot [0.08:2][1.4:1.9] \
"pbp-b06.dat" i 0 u 4:5:6 w yerrorlines lc "violet" pt  4 ps 1 title columnhead, \
"pbp-b06.dat" i 1 u 4:5:6 w yerrorlines lc "violet" pt 12 ps 1 title columnhead, \
"pbp-b06.dat" i 2 u 4:5:6 w yerrorlines lc "violet" pt  6 ps 1 title columnhead

quit

