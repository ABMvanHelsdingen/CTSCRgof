# CTSCRgof

This package contains functions to implement all GOF methodologies presented in "Goodness of Fit tests for Continuous-Time Spatial Capture-Recapture Models" (2026\*).<br/>

The frequentist time rescaling tests can be computed with *timeRescaling()* and the residuals with *getResiduals()*. Both functions are in *Frequentist.R*
A generic function for Bayesian Posterior Predictive Checks (PPCs) can be found in *Bayesian.R*. A metric must be specified; those used in the Appendix can be found in *Metrics.R*.<br/>

A function for the hazard function of the continuous-time SCR model must be provided for all GOF tests. *Hazards.R* includes hazard functions for three SCR models mentioned in the text. The Bayesian tests also require a simulation function. The file *Simulation.R* includes simulation functions for the same three models.<br/>


The vignette simulates a SCR model and demonstrates all the frequentist GOF tests.

