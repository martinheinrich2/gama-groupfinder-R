# gama-groupfinder

This code was used in an experiment to test if there is an empirical natural linking
lenght for linking groups of galaxies.
A simple percolation or friends-of-friends (FOF) algorithm was split into two tasks.
First link galaxy pairs within a given line-of-sight and transverse distance, then
determine groups of galaxies from these pairs.
The FOF algorithm usually is calibrated with mock data to produce galaxy group catalogues.
Since there is no overall consensus about the optimal linking lenghts in the literature,
this code can be run for many different linking-lengths. The distribution of groups, their
richness varies in the process. It is possible that a natural empirical linking lenght
would show some stability with variations for small changes in the linking length.

Two data sets were used, the science targets from the GAMA-Survey [TilingCat input catalogue](http://www.gama-survey.org/dr3/data/cat/EqInputCat/v46/TilingCat.fits) and redshifts from [DistanceFrames](http://www.gama-survey.org/dr3/data/cat/LocalFlowCorrection/v14/DistancesFrames.fits).
They were combined and a volume-limited sample was created using [TOPCAT](http://www.star.bristol.ac.uk/~mbt/topcat/).

## Procedure

1. Calculates line-of-sight (i.e. radial) comoving distance in Mpc, transverse comoving distance in Mpc,
   line of velocity in km/s and zeta. Converts RA from degree to radians and DEC from degree to radians. Writes 
   selected columns to file: "reduced_converted_sample.csv" 
   and full data to: "converted_sample.csv"

2. Galaxy-pair-finder.R finds galaxy pairs within a given range
   of line-of-sight (in km/s) and transverse (in Mpc) values. It writes results as linknumbner, CATAID1, CATAID2.
   It uses data.table instead of data.frame to speed the code up.

3. Galaxy-group-finder-final.R finds galaxy groups from galaxy pairs. Loops over a folder with
   galaxy-pairs, links them to groups, writes groups and richness (distribution of groups) to csv-files.

4. The poweRlaw_all_tail_final.R calculates powerlaw-fit and parameters. Loops over all files in
   folder richness and writes slope (alpha), sigma (standard error), x_min (data value beyond which
   the distribution was fitted), normalisation (C), tailRatio (tail length/all elements).
   It uses the poweRlaw package.

5. Analyse powerlaw parameters poweRlaw_plot.R. Creates heatmaps from the poweRlaw 
   results file: los, trans and tail-ratio, variable tail-ratio, alpha, normalisation, xmin, maxrichness
