
set term postscript eps enhanced color
set output "Cs3.eps"
set encoding iso_8859_1

unset surface
set nogrid
set view map
set contour base
set cntrlabel format '%8.3g' font ',7' start -0.001 interval 0.001
set cntrparam levels incremental -0.1,0.002,0.01
unset clabel
set dgrid3d 100,100,6
set xrange[6.7:19]
set yrange[6.7:19]
set pm3d
set ztics 0.01

# Define cores e paleta
set linetype 1 lc rgb "white"  # Linhas de contorno brancas
set palette defined (0 "#440154", 0.2 "#31688e", 0.4 "#35b779", 0.6 "#fde725", 1 "#ff0000")

# Configurar a barra de cores (colorbar)
set colorbox  # Garante que a barra de cores está visível
set cbrange [-0.05:0.2]  # Intervalo explícito para a barra de cores
set cblabel "{/NewCenturySchlbk-Bold=12 E/E_{h}}" offset 0,-0.5  # Unidade na barra de cores
# Adicionar o valor mínimo como um rótulo personalizado na barra de cores
#set cbtics add (-0.0222 "Min: -0.0222 E_h")  # Etiqueta do mínimo

# Rótulos dos eixos
set xlabel "{/NewCenturySchlbk-Bold=12 r_{12}/{a_{0}}"
set ylabel "{/NewCenturySchlbk-Bold=12 r_{13}/{a_{0}}}"
set zlabel "{/NewCenturySchlbk-Bold=12 E/E_{h}}"

# Plotar superfície e contornos
splot 'Cs3eq.dat' using 1:2:4 with pm3d, \
      '' using 1:2:4 with lines lt 1 lc rgb "white" notitle
