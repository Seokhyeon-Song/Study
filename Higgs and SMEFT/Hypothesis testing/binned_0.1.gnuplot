cd 'C:\Users\hypercube256\Desktop\Study\Higgs and SMEFT\Hypothesis testing'
set terminal pdfcairo font "Times-New-Roman,12" size 4,3
set output "binned_0.1GeV.pdf"
set datafile separator ","
set dgrid3d 200,200 gauss 0.0000265,0.0001
set table 'test.dat'
splot 'fitresult_binned_0.100000.csv'
unset table

set contour base
set cntrparam levels discrete 1000,3000,10000,30000,100000,300000
unset surface
set table 'cont.dat'
set format cb "%.s%c"
set format z "%.s%c"
splot 'fitresult_binned_0.100000.csv'
unset table

reset
unset key
set palette rgbformulae -33,-13,-10
set logscale zcb
set cbrange [100:10000000]
set xrange [0:0.0531]
set yrange [0.4:0.6]
set xlabel "Coupling constant"
set ylabel "Loop particle mass (GeV)"
set format cb "%.s%c"
plot 'test.dat' with image, 'cont.dat' every ::40:1::1 with lines lt -1 lw 1.5, 'cont.dat' every :::1:20:1 with lines lt -1 lw 1.5, 'cont.dat' every :::::0 with lines lt -1 lw 1.5, 'cont.dat' every :::2 with lines lt -1 lw 1.5, 'cont.dat' every ::30:1:30:1 with labels
