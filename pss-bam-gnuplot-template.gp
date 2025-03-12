set terminal postscript eps color "Arial" 14;
#############################################################
#                CUSTOMIZE THESE OPTIONS                    #
#############################################################
# This template makes a map-damage like plot from the
# output file of pss-bam.pl run from merged reads, i.e.,
# using the -m option to pss-bam.pl
# Customize these 4 lines to control what the input file and
# output files are and how much damage was present in the
# reads. Then run:
# gnuplot THIS_FILE
# to generate the plot
#############################################################
set ytic 0.05;                        # tic marks on y-axis #
set output "plot.eps";          # output filename #
datafile = "filename.pss.counts.txt"; # input .pss.counts.txt file #
max_damage = 0.1           # high fraction of C->T in plot  #
#############################################################

set multiplot;
set size 0.5, 1;
set grid;
set style line 1 lt 1 lc 7 lw 0.5;
set style line 2 lt 1 lc rgb "#8b0000" lw 3;
set style line 3 lt 1 lc rgb "#228b22" lw 3;
set style line 4 lc rgb "green"; # A = green in Sanger
set style line 5 lc rgb "blue";  # C
set style line 6 lc rgb "black"; # G
set style line 7 lc rgb "red"; # T
set nokey;
unset xtic;
set tmargin 0.9;
set bmargin 0.5;
set style data histogram;
set style histogram rowstacked;
set style fill transparent solid 0.4 noborder;
set boxwidth 0.8;
set style rect fc lt -1 fs solid 0.35 noborder;
set obj rect from -1.0, graph 0 to 1.5, graph 1;
r = 14;
set origin 0, 0;

#set title "plot-title";

plot [-1.0:17.0][0:max_damage] datafile index 0 using ($2+$6+$10+$14) axes x1y2 ls 4,\
     '' index 0 using ($4+$8+$12+$16) axes x1y2 ls 6,\
     '' index 0 using ($3+$7+$11+$15) axes x1y2 ls 5,\
     '' index 0 using ($5+$9+$13+$17) axes x1y2 ls 7,\
     '' index 0 u ($1+2):($3==0?NaN:$3/($3+$7+$11+$15)) t "C->A" w l ls 1,\
     '' index 0 u ($1+2):($4==0?NaN:$4/($4+$8+$12+$16)) t "G->A" w l ls 3,\
     '' index 0 u ($1+2):($5==0?NaN:$5/($5+$9+$13+$17)) t "T->A" w l ls 1,\
     '' index 0 u ($1+2):($6==0?NaN:$6/($2+$6+$10+$14)) t "A->C" w l ls 1,\
     '' index 0 u ($1+2):($8==0?NaN:$8/($4+$8+$12+$16)) t "G->C" w l ls 1,\
     '' index 0 u ($1+2):($9==0?NaN:$9/($5+$9+$13+$17)) t "T->C" w l ls 1,\
     '' index 0 u ($1+2):($10==0?NaN:$10/($2+$6+$10+$14)) t "A->G" w l ls 1,\
     '' index 0 u ($1+2):($11==0?NaN:$11/($3+$7+$11+$15)) t "C->G" w l ls 1,\
     '' index 0 u ($1+2):($13==0?NaN:$13/($5+$9+$13+$17)) t "T->G" w l ls 1,\
     '' index 0 u ($1+2):($14==0?NaN:$14/($2+$6+$10+$14)) t "A->T" w l ls 1,\
     '' index 0 u ($1+2):($15==0?NaN:$15/($3+$7+$11+$15)) t "C->T" w l ls 2,\
     '' index 0 u ($1+2):($16==0?NaN:$16/($4+$8+$12+$16)) t "G->T" w l ls 1;

set origin 0.5, 0;
unset obj;
set obj rect from 14.5, graph 0 to 17, graph 1;
plot [-1.0:17.0][0:max_damage] datafile index 1 using ($2+$6+$10+$14) axes x1y2 ls 4,\
     '' index 1 using ($4+$8+$12+$16) axes x1y2 ls 6,\
     '' index 1 using ($3+$7+$11+$15) axes x1y2 ls 5,\
     '' index 1 using ($5+$9+$13+$17) axes x1y2 ls 7,\
     '' index 1 u ($1+2):($3==0?NaN:$3/($3+$7+$11+$15)) t "C->A" w l ls 1,\
     '' index 1 u ($1+2):($4==0?NaN:$4/($4+$8+$12+$16)) t "G->A" w l ls 3,\
     '' index 1 u ($1+2):($5==0?NaN:$5/($5+$9+$13+$17)) t "T->A" w l ls 1,\
     '' index 1 u ($1+2):($6==0?NaN:$6/($2+$6+$10+$14)) t "A->C" w l ls 1,\
     '' index 1 u ($1+2):($8==0?NaN:$8/($4+$8+$12+$16)) t "G->C" w l ls 1,\
     '' index 1 u ($1+2):($9==0?NaN:$9/($5+$9+$13+$17)) t "T->C" w l ls 1,\
     '' index 1 u ($1+2):($10==0?NaN:$10/($2+$6+$10+$14)) t "A->G" w l ls 1,\
     '' index 1 u ($1+2):($11==0?NaN:$11/($3+$7+$11+$15)) t "C->G" w l ls 1,\
     '' index 1 u ($1+2):($13==0?NaN:$13/($5+$9+$13+$17)) t "T->G" w l ls 1,\
     '' index 1 u ($1+2):($14==0?NaN:$14/($2+$6+$10+$14)) t "A->T" w l ls 1,\
     '' index 1 u ($1+2):($15==0?NaN:$15/($3+$7+$11+$15)) t "C->T" w l ls 2,\
     '' index 1 u ($1+2):($16==0?NaN:$16/($4+$8+$12+$16)) t "G->T" w l ls 1;
