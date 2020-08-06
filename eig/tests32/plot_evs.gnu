
set terminal pdfcairo size 15in, 9in

# set aspect ratio to 1
set size ratio -1

# set output file name
set output "L32beta05.pdf"

# plot to sets of eigenvalues together on the same plot
plot [-0.5:2.5] [-1:1] "b05000k00000-HF.eig", "b05000k00000-OV.eig"

quit

