# EDPcausal

The codes are to run two simulation scenarios in the main manuscript. 

# Explain the simulation scenarios

We consider two different simulation scenarios. In scenario 1, the sets of covariates for
both confounders and effect modifiers are identical (X = U). In scenario 2, the set of effect
modifiers includes more covariates than the set of confounders (i.e.,X   U; a proper subset).
Table 1 describes the data generating mechanism for scenario 1. For scenario 2, we consider
two additional effect modifiers X6 and X7.

# MCMC running

    - Run 'sim_new_s1.R' or 'sim_new_s2.R' script for obtaining posterior samples for scenario 1 or scenario 2, respectively.
    - Executing 'colMeans(temp)' produces the cluster-specific effects.
    
