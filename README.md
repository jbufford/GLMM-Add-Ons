# GLMMs
These scripts provide add-on code for running GLMMs in R.  

<b>AIC Tables.R</b>: Runs a set of <i>a priori</i> models and makes a delta AIC or delta AICc table from the results.  It also returns the suite of models.  It works well with most models, but the AIC values may not be correct with some correlation structures in nlme (e.g. varPower).

<b>Model Checks JB.R</b>: These run a suite of model checks on a given (generalized) linear mixed effects model.  These include graphs of raw data, residual plots, tests for overdispersion, catterpillar plots, etc.  Optionally runs tests of influence using code modified from the package influence.ME to allow for more complex model structures and to run using parallel processing (<b>InfluenceJB.R</b>).  This code also calls <b>rsquaredglmm</b> (https://github.com/jslefche/rsquared.glmm).  To get r-squared values, the user must also access the rsquaredglmm code.  The whole model check output can be exported to a .pdf using R Markdown (<b>Model_Checks_JB.Rmd</b>).  Figures only can be exported to a pdf with a pdf device.

<b>Graph GLMMs.R</b>: Code for bootstrapping, calculationg and graphing the coefficients and 95% CI from a (g)lmm.  This script includes a function for use in the function boot (package boot) for bootstrapping fixed effect coefficients from lmer or nlme models.  The plot.coef function takes as its input either a data.frame of MCMC CI or a list of boot.ci-generated objects, or a boot object from boot().  It generates ggplots of coefficients and 95% CI, with variables that don't overlap 0 indicated as filled circles.  Theme properties can be modified in the final ggplots.  The input data structure for MCMC CI is currently inflexible.

<b>GLMM Simulations.R</b>: Runs simulations to test the effect of heteroscedastic variances on coefficient estimates.  This code was purpose-made for a specific test and is still incomplete.  It has not been tested on a broader variety of possible models.
