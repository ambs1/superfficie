reset
#set nokey

set term postscript eps enhanced color
set output "Cs3.eps"
set encoding iso_8859_1

unset surface
set nogrid
set view map
set contour base
set cntrlabel format '%8.3g' font ',9' start -0.001 interval 0.001
set cntrparam levels incremental -0.1,0.002, 0.01
unset clabel
set dgrid3d 100,100,6
#set xrange[6:20]
#set yrange[6:20]

set ztics 0.01
set xlabel "{/NewCenturySchlbk-Bold=12 r_{12}/{a0}}"
set ylabel "{/NewCenturySchlbk-Bold=12 r_{13}/{a0}}"
set zlabel "{/NewCenturySchlbk-Bold=12 E/E_{h}}"

# Altere a cor aqui (ex: "black", "red", "dark-green", etc)
splot 'Cs3eq.dat' using 1:2:4 title "" with lines lc rgb "black"
pause -1 "Fecha?"
