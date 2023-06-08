The code in here calculates the fluctuations in the UE from some input file containing events.

Here are the instructions for using it:

Open ROOT and load getFluctuations_sectorExclude.C
    .L getFluctuations_sectorExclude.C
run getFluctuations over your root file containing the events.
    getFluctuations(<# of window sizes>, <array of window sizes>, <input file>, <output file>, %<subtraction level>, %<allow partial exclusions>, %<# of events>, %<Et>, %%<towers to exclude>)
Agruments marked with a % are optional. Arguments marked with %% are recommended not to be changed. More details on the function's arguments are given in the function definition.

This produces a single root file named <ouput filet>_unsubtracted.root or <output file>_subtracted_noflow.root or <output file>_subtracted_withflow.root
This file contains histograms of the STD of window energy in an event vs. the event centrality for each inputted window size for each of the EMCal, IHCal, OHCal, and total calorimeter.
Additionally, the file contains histograms of STD window energy vs. total event energy and event impact parameter, and the average evergy in a window vs. these quantities.

If you wish to save the histograms produced in the first step as pdfs, load and run fluctuation_plotter_strangesizes.C, using the output from getFluctuations() as an input.
    .L fluctuation_plotter_strangesizes.C
    fluctuation_plotter_strangesizes(<# of window sizes>, <array of window sizes>, <input file>, <output directory>)
    
Finally, load and run STDbargraphs_nofit.C using the root rile produced by getFluctuations() as the <input file>.
    .L STDbargraphs_nofit.C
    STDbargraphs(<input file>, <output directory>)

This function produces the final plots of the STD window energy averaged over all events as a function of window size, with events grouped into centrality ranges.
The final plots are printed as pdfs into <output directory>.
As an intermediate step, the funciton projects the histograms produced by getFluctuations() onto the y-axis, producing 1D histograms of STD energy for each window size and calorimeter.
These intermediate histograms are also printed as pdfs.
    

The process for running producing these plots may be made more straightforward in the future.
