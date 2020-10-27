# CausalAnalysisforSemiCompRisks
This is an R-package for analyzing semi-competing risks data using causal mediation inference. This package is aimed to illustrate the methodology presented in "Semiparametric causal mediation modeling of hepatitis on mortality through liver cancer incidence".

## Installation
```r
devtools::install_github("eric40065/CausalAnalysisforSemiCompRisks")
```

## Usage
The following code can reproduce everything mentioned in our paper. It will take one hour to reproduce Figure 1 and 2 and over 10 days to reproduce Table 1.
```r
library(CausalAnalysisforSemiCompRisks)
# To reproduce Figure 1(a), 1(b), 1(e), 1(f), 1(i), and 1(j) in page 19
simulation(1, hypo = 'null')
# To reproduce Figure 1(c), 1(d), 1(g), 1(h), 1(k), and 1(l) in page 19
simulation(1, hypo = 'alter')
# To reproduce the first part of Table 1 in page 21
Table11 = simulation(2, hypo = 'null', sample_size = 300)
# To reproduce the second part of Table 1 in page 21
Table12 = simulation(2, hypo = 'null', sample_size = 1000)
# To reproduce the third part of Table 1 in page 21
Table13 = simulation(2, hypo = 'alter', sample_size = 300)
# To reproduce the last part of Table 1 in page 21
Table14 = simulation(2, hypo = 'alter', sample_size = 1000)
# To reproduce Figure 2(a)--(d) in page 25
HBV_result = CASCR(REVEAL_HBV, get_variance = c('asymptotic', 'bootstrap'), plot_result = T)
# To reproduce Figure 1(a)--(j) in page 20 of supplement material
HBV_sen_ana = CASCR(REVEAL_HBV, sen_ana = T, get_variance = NULL)
```
The following code provides a glimpse of the analysis.

### simulation
```r
library(CausalAnalysisforSemiCompRisks)
# The result of the unbiasedness as presented in Figure 1 in our paper. We repeat it 1,000
# times to get more accurate result.
simulation(1, hypo = 'null', repeat_size = 100)
simulation(1, hypo = 'alter', repeat_size = 100)

# The result of coverage rate as presented in Table 1 in our paper. We repeat it 1,000 times
# and we set get_variance as c('aymptotic', 'bootstrap') to get more accurate result.
# This, however, spend plenty of time.
simulation(2, hypo = 'null', sample_size = 1000, repeat_size = 100, get_variance = 'asymptotic')
simulation(2, hypo = 'alter', sample_size = 1000, repeat_size = 100, get_variance = 'asymptotic')
```

### REVEAL-HBV
#### Abstract
- REVEAL-HBV is a community-based prospective cohort study conducted in Taiwan.
- We adjusted age at the cohort entry and the history of alcohol consumption (yes vs. no) as covariates in our model.
#### Dictionary
REVEAL-HBV contains 11,946 male patients with 7 columns.
- hcc.time: time to liver cancer incidence or censored time (days).
- die.time: time to death or censored time (days).
- hcc.case: whether or not the liver cancer incidence is observed (0: no vs. 1:yes).
- dieall.case: whether or not death is observed (no: 0 vs. yes: 1).
- HBSAG: hepatitis B surface antigen (negative: 0 vs. positive : 1)
- AGE: age at the cohort entry
- alcohol1: the alcohol consumption history (no: 0 vs. yes: 1)
#### Code
```r
library(CausalAnalysisforSemiCompRisks)
# The result of analyzing the REVEAL-HBV dataset as discussed in Section 7 in our paper.
# We set downsample as 1 and get_variance as c('aymptotic', 'bootstrap') in our paper.
result = CASCR(REVEAL_HBV, downsample = 5, get_variance = 'asymptotic', plot_result = T)

# To get a quick result including bootstrap variance, one can run the following code.
# It will take around 15 minutes.
result = CASCR(REVEAL_HBV, downsample = 60, get_variance = c('asymptotic', 'bootstrap'), plot_result = T)

# The result of sensitivity analysis as presented in Section 11 of Supplement Material.
result = CASCR(REVEAL_HBV, sen_ana = T, get_variance = NULL)
```

## More details
More details such as the input data type, parallel programming, and the output can be found using
```r
?CASCR
```
