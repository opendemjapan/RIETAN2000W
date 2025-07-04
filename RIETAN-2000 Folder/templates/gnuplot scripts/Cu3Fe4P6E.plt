# Script file for observed, calculated, and difference patterns

# Destination of the output \
iout = 0: Screen (plot only observed and calculated intensities). \
iout = 1: Screen (complete graph). \
iout = 2: Encapsulated PostScript (EPS) file.
iout = 1

if (iout <= 1) \
set terminal windows color "Helvetica" 14; \
set size; \
set output

if (iout == 2) \
   set terminal postscript eps enhanced monochrome solid "Times-Roman" 26; \
   set size 2.0, 1.25; \
   set output 'Cu3Fe4P6E.eps'

xmin = 10 # Minimum for the x-axis
xmax = 110 # Maximum for the x-axis
ymin = -5000 # Minimum for the y-axis
ymax = 45000 # Maximum for the y-axis
offset_delta =-2000 # Offset of the difference curve along the y axis
len_bar = 0.010 * (ymax - ymin) # Length of bars for peak positions
y_phase1 = 4000 # y coordinate of bars for phase No. 1
y_phase2 = 2700 # y coordinate of bars for phase No. 2
y_phase3 = 1400 # y coordinate of bars for phase No. 3
set ytics 0, 10000, 40000 # Start, increment, and end values of tics for the  y axis
set pointsize 0.7 # Size of '+' symbols for observed intensities
set mxtics 2 # minor tics (midpoints) for the x-axis
set mytics 2 # minor tics (midpoints) for the y-axis

if (iout <= 1) \
   set xlabel '2-theta/deg' 0.0, 0.2 'Helvetica, 16'; \
   set ylabel 'Intensity' 0.5, 0.0  'Helvetica, 16'

if (iout == 2) \
   set xlabel '2{/Symbol q} / deg' 0, 0.20 'Times-Roman,30'; \
   set ylabel 'Intensity' 0.25, 0 'Times-Roman,30'

# To learn colors and line types, type 'test' in the console window
if (iout == 0) plot [*:*] [*:*] \
   'Cu3Fe4P6E'using 1:3 notitle with lines linetype 11 linewidth 1, \
   '' using 1:2 notitle with points 2

if (iout == 1) plot [xmin:xmax] [ymin:ymax] \
   'Cu3Fe4P6E.pat' using 1:3 notitle with lines linetype 11 linewidth 1, \
   '' using 1:2 notitle with points 2, \
   '' using  1:($2 - $3 + offset_delta) notitle with lines linetype 4 linewidth 1, \
   '' using  8:(y_phase1):(len_bar) notitle with yerrorbars linetype 3 linewidth 1 pointtype 0
   
if (iout == 2) plot [xmin:xmax] [ymin:ymax] \
   'Cu3Fe4P6E.pat' using 1:3 notitle with lines linetype 5 linewidth 1, \
   '' using 1:2 notitle with points 1, \
   '' using  1:($2 - $3 + offset_delta) notitle with lines linetype 3 linewidth 1, \
   '' using  8:(y_phase1):(len_bar) notitle with yerrorbars linetype 2 linewidth 1 pointtype 0, \
   '' using 14:(y_phase2):(len_bar) notitle with yerrorbars linetype 4 linewidth 1 pointtype 0, \
   '' using 20:(y_phase3):(len_bar) notitle with yerrorbars linetype 4 linewidth 1 pointtype 0 
