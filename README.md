# Bayesian estimation of the GJR(1,1) model with student-t innovations and explanatory variables
The GJR time series model is a proficient model for time series analysis. Furthermore, Bayesian estimations of these model is essential for small data sets analysis. Here, we provide a R code for estimating the GJR(1,1).
The particularity of this code consists in the introduction of explanatory variables.

First of all, the file **Model, Method and Application.pdf** provides three elements: 1) A clear definition of the GJR(1,1) model with explanatory variables 2) A gentle introduction to the Metropolis sampler 3) A direct application on financial time series of the Bayesian estimate of the GJR(1,1) model with explanatory variables.
Secondly, the folder **code** contains the R file which gives an R API to the GJR models and an examples of use. In the subfolder **models**, the curious reader may take a look at the stan files coding the GJR models and the Metropolis-Hasting algorithm. In addition, the subfolder **application** contains the R files and the financial data used for the application of the GJR(1,1) model with explanatory variables.
