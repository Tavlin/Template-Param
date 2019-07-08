# Template-parametrisation
Code for my bachelor thesis.
Fitting and plotting of MC templates on data

1. ## Prerequisite
  * ### Normal afterburner output
     This is simply the normal output of the afterburner which we used for RUN 1
     and 2. Basically what you would get, when you are doing a normal analysis
     with the afterburner.

     This is the main output file which will be used in this analysis.

  * ### Lower & higher rebinning
     We normally rebin the invariant mass-axis, and this rebinning has to be
     changed into a lower and a higher rebinning separately. This is one of the
     systematic variations and is therefore needed for the calculation of the
     systematic uncertainties.

     This change needs to be done by hand in the afterburner to obtain the
     needed output!

  * ### No Rebinng
     No rebinning (or the "1-rebinning" how I call it) means you change the
     value in the invariant mass-axis rebinning to 1, so it does not get
     rebinned. This is needed for the combination of the correlated background
     templates, which may be needed because of the statistical uncertainties.

     This change needs to be done by hand in the afterburner to obtain the
     needed output!

  * ### Smaller and wider background fitting
     For systematic uncertainties this variation is normally included in the
     afterburner however, the afterburner changes the fitting range from right
     of the peak to left of the peak for pi0. Since we use this template method
     to be able to describe the pions where at least one of the two photons
     converted into e+e- outside of our PID. This part of the pi0 peak should be
     in the same mass-range where the background fitting range variation of the
     afterburner is.
     So I decided to change the systematic variation into broadening and
     narrowing the fit range of the background instead of translating it to
     other invariant mass-values.

     This change needs to be done by hand in the afterburner to obtain the
     needed output!
