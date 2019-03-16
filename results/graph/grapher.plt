set xtics format
set logscale x
set ytics format
set logscale y
set grid y my
set grid x mx
set xlabel "Periods (T)" font "Arial, 6"
set ylabel "PSA (m/s2)" font "Arial, 6"
set title "Comparative of Atennuations Laws"
plot "toroE02spectrum.dat" using 1:2 title "toroE2002" with lines lc "green" lw 2,\
"toroM02spectrum.dat" using 1:2 title "toroM2002" with lines lc "purple" lw 2,\
"ab2006spectrum.dat" using 1:2 title "ab2006" with lines lc "red" lw 2,\
"ab2011spectrum.dat" using 1:2 title "ab2011" with lines lc "blue" lw 2,\
"pezeshk2011spectrum.dat" using 1:2 title "pezeshk2011" with lines lc "yellow" lw 2,\
"silva2002.dat" using 1:2 title "silva2002" with lines lc "black" lw 2,\
"campbell2003spectrum.dat" using 1:2 title "campbell2003" with lines lc "gray" lw 2