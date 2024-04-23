# EDPcausal

The codes are to run two simulation scenarios in the main manuscript "*Bayesian nonparametric model for heterogeneous treatment effects with zero-inflated data*" by Chanmin Kim, Yisheng Li, Ting Xu, and Zhongxing Liao.  

# Explain the simulation scenarios

We examine two distinct simulation scenarios. In Scenario 1, the degree of violation of the overlap assumption is moderate, while in Scenario 2, the violation of the overlap assumption is severe. In Scenario 1, there is overlap between the two groups when considering the values of X2 within the range of -1 to 2. However, the overlap region becomes significantly narrower in Scenario 2. For the simulated data (n = 600), we generate one binary covariate, one discrete covariate, and three continuous covariates as effect modifiers. The binary treatment is generated using a logistic regression model, while the outcome is generated from a truncated normal distribution. The probability π for generating excessive 0 values is derived from the logistic regression model.

# MCMC running

* Scenario 1: 
  1. Download the `sim1` folder. 
  2. run `sim_edp.R` for the proposed model; run `sim_bcf.R` or `sim_dpglm.R` for the competing models.
* Scenario 1: 
  1. Download the `sim2` folder. 
  2. run `sim_edp.R` for the proposed model; run `sim_bcf.R` or `sim_dpglm.R` for the competing models.
    
