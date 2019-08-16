# Supervised Randomization

Code accompanying the paper: 
Haupt, Jacob, Gubela & Lessmann (2019). "Affordable Uplift: Supervised Randomization in Controlled Experiements".

## Summary
Supervised randomization integrates the targeting model into randomized controlled trials and A/B tests as a stochastic policy. 
By replacing random targeting with supervised randomization, we reduce the costs of running randomized experiments or allow continuous collection of randomized trial data. We show how to fully correct downstream analysis for the bias effect by the supervised treatment allocation with the true probabilites to receive treatment, which are logged at the time of targeting. 

## Data
The data we use for empirical monte carlo simulation is public at https://archive.ics.uci.edu/ml/datasets/Bank%2BMarketing