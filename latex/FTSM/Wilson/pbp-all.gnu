
set terminal pdfcairo size 4in, 3in

# set aspect ratio to 1
#set size ratio -1

# set output file name
set output "pbp-all.pdf"

#set title "Wilson fermions, {/Symbol b} = 2.0 (blue) and {/Symbol b} = 6.0 (violet)"
#set xlabel "T / {/Symbol m}"
set xlabel "T = 1 / n_t"
set ylabel "<~{/Symbol Y}{.8-}{/Symbol Y}>"

set logscale x

# pt = pointtype (12 = diamonds); ps = pointsize
plot [0.05:1][1.4:1.9] \
"pbp-b02.dat" i 0 u 3:5:6 w yerrorlines lc "blue" pt  4 ps 1 title columnhead, \
"pbp-b02.dat" i 1 u 3:5:6 w yerrorlines lc "blue" pt 12 ps 1 title columnhead, \
"pbp-b02.dat" i 2 u 3:5:6 w yerrorlines lc "blue" pt  6 ps 1 title columnhead, \
"pbp-b06.dat" i 0 u 3:5:6 w yerrorlines lc "violet" pt  4 ps 1 title columnhead, \
"pbp-b06.dat" i 1 u 3:5:6 w yerrorlines lc "violet" pt 12 ps 1 title columnhead, \
"pbp-b06.dat" i 2 u 3:5:6 w yerrorlines lc "violet" pt  6 ps 1 title columnhead
#"pbp-b08.dat" i 0 u 3:5:6 w yerrorlines lc "red" pt  4 ps 1 title columnhead

quit

