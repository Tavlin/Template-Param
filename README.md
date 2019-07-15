Template-parametrisation
Code for my bachelor thesis.
Fitting and plotting of MC templates on data

1. ## Prerequisite
  * ### Normal afterburner output
     This is simply the normal output of the afterburner which we used for RUN 1
     and 2. Basically what you would get, when you are doing a normal analysis
     with the afterburner.

     This is the main output file which will be used in this analysis.

  * ### Lower & higher rebinning
     We usually rebin the invariant mass-axis. This rebinning has to be
     done separately into a lower and a higher rebinning separately. This is one of the
     systematic variations and is therefore needed for the calculation of the
     systematic uncertainties.

     This Setting is needed to be changed by hand!

  * ### No Rebinng
     No rebinning (or the "1-rebinning" how I call it) means you change the
     rebinning value in the invariant mass-axis rebinning to 1, so it does not
     get rebinned. This is needed for the combination of the correlated
     background templates, which is needed since the correlated background
     templates are used over a larger _p_ T range, then one _p_ T bin. So to be
     able to combine them they need to be binned in the same way.

     This Setting is needed to be changed by hand!

  * ### Smaller and wider background fitting
     For systematic uncertainties this variation is usually included in the
     afterburner however, the afterburner changes the fitting range from the
     right to the left of the peak. Since we use this template method to
     describe the pions where at least one of the two photons converted into
     e+e- outside of our PID, this part of the pi0 peak should be in the same
     mass region as the fitting range variation of the afterburner.
     So I decided to change the systematic variation into broadening and
     narrowing the fitting of the background instead of translating it to other
     invariant mass-values.

     This Setting is needed to be changed by hand!

  * ### Setting up a Files.txt
     To run the code smoothly I wrote a bash script that reads all the needed
     path variables (like where the different output files are). Those should be
     contained in a single file, called _Files.txt_
     In the bash script are the different needed variables listed as follows:

     > $a == what should be plotted

     > $b == SaveFile Format (i.e. eps, pdf, png etc. without a dot!!)

     > $c == path to the ESD MC

     > $d == path to the ESD data

     > $e == path to the normal framework output

     > $f == path to the framework output without rebinning

     > $g == path to the framework output with higher rebinning

     > $h == path to the framework output with lower rebinning

     > $i == path to the framework output with narrow uncorr back fit range

     > $j == path to the framework output with wide uncorr back fit range

     > $k == cut string (since the cut string log can contain more then one and this is made only for one cut string!)

     > $l == SaveAppendix (i.e. pp13TeV, or pp2016)


  2. ## Running the analysis
     To run the analysis after you set up everything simply use the command:
     ```bash
     $ bash templateanalysis.sh
     ```

     By default the bash script runs the full analysis, with systematic
     uncertainty calculation and plotting. If you wish to only do some parts,
     you will have to change a few things in the bash script or comment some
     lines out by hand.

     By default the folders in which the plots are saved will be deleted and new
     ones created, keep that in mind!
