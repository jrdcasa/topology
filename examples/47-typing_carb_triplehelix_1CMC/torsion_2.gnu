reset
###############################################################################
set term wxt 1 enhanced dashed size 500,400 font "Arial,10"
set multiplot layout 1,1
set encoding iso_8859_1

set style line 1 lt 1 ps 0.4 lc rgb "black"  pt 4 lw 2.0
set style line 2 lt 2 ps 0.4 lc rgb "blue"  pt 4 lw 2.0
set style line 3 lt 2 ps 0.4 lc rgb "blue"   pt 4 lw 2.0
set style line 4 lt 1 ps 0.4 lc rgb "green"  pt 4 lw 2.0
set style line 5 lt 2 ps 0.4 lc rgb "yellow" pt 4 lw 2.0
set style line 6 lt 2 ps 0.4 lc rgb "orange" pt 4 lw 2.0

#set xrange[-pi:pi]
set xrange[-360:360]
set xlabel "{/Symbol F} (rad)"
set ylabel "E_{tors} (kcal/mol)"
set key top center

# Proper multiple
set angle degrees
t1(x)  = 2.677760*(1+cos(1*x-180)) + 0.125520*(1+cos(2*x-180)) + 2.552240*(1+cos(3*x-0))
t2(x)  = 2.284880*(1+cos(1*x-0)) + 1.213360*(1+cos(2*x-0)) + 1.799120*(1+cos(3*x-0))

p t1(x)  title 'T1' w l ls 1,\
  t2(x)  title 'T2' w l ls 2


unset multiplot

