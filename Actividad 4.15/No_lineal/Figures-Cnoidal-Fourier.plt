#
# Plot results using Gnuplot
#

reset

BaseOutput = 0
Output = BaseOutput

LW = 1
if (Output > 0) LW = 3.5
set border lw 1.4

Xscale = 1; Yscale = 0.5
set style line 1 linetype 1 linecolor rgb "blue" lw LW
set style line 2 linetype 1 linecolor rgb "red" lw LW
set style line 3 linetype 1 linecolor rgb "blue" lw LW dt(4,4)
set style line 4 linetype 1 linecolor rgb "red" lw LW dt(4,4)
set style line 6 linetype 6 linecolor rgb "black" ps 0.6 lw 2.5

Dir_Fourier = "C:/JF/Internet/Steady-waves/"
Dir_Cnoidal = Dir_Fourier."Cnoidal/"

#####################################################################
# Generic Surface profile - testing current versions
#####################################################################

set key bottom center reverse Left horizontal height 2 width +6
set format x "%3.0f"
set format y "%4.1f"
if (Output==0) {set xlabel "x/d"; set ylabel "y/d"}
if (Output>0) {
    set xlabel "$x/d$" offset 0,-0.5; 
    set ylabel "$\\dfrac{\\eta}{d}$" offset -1
    }
set xtics nomirror
set ytics 0, 0.5, 2
set yrange [0:1.6]
set autoscale x

title = "Cnoidal and Fourier"

File=Dir_Cnoidal."Surface-Test"

load Dir_Fourier.'SetOutput.plt'
if (Output >= 0) plot Dir_Fourier."Surface.res" using 1:2 title "Fourier approximation" with lines ls 1,\
     Dir_Cnoidal."Surface.res" using 1:2 title "Cnoidal theory" with lines ls 4;\
     pause -1 "Hello Surface - testing current versions"

############################################
# Velocity profiles
############################################

Xscale = 1; Yscale = 0.6

Output = BaseOutput 

set key top center reverse Left horizontal width +4
set xtics offset 0,-0.5
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) {set ylabel "y/d"}
if (Output>0) {
    set ylabel "$\\dfrac{y}{d}$" offset -1
    }
if (Output==0) set xlabel "Velocity (dimensionless w.r.t. g \\& d)";
if (Output>0) set xlabel "$(u,v)/\\sqrt{g d}$" offset 0,-1;
set title "Cnoidal and Fourier: Velocity profiles over half a wave"
# set title ""
File = Dir_Cnoidal.'uv'
load Dir_Fourier.'SetOutput.plt'
if (Output >= 0) plot Dir_Fourier."Flowfield.res" using 2:1 title "$u$ - Fourier" with lines ls 1,\
     Dir_Cnoidal."Flowfield.res" using 2:1 title "$u$ - Cnoidal" with lines ls 3,\
     Dir_Fourier."Flowfield.res" using 3:1 title "$v$ - Fourier" with lines ls 2,\
     Dir_Cnoidal."Flowfield.res" using 3:1 title "$v$ - Cnoidal" with lines ls 4;\
pause -1 "Hello uv"

unset output
