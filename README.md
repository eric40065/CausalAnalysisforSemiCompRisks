# CausalAnalysisforSemiCompRisks
This is an R-package for analyzing semi-competing risks data using causal mediation inference. This package is aimed to illustrate the methodology presented in "Semiparametric causal mediation modeling of hepatitis on mortality through liver cancer incidence".

## Installation
```r
devtools::install_github("eric40065/CausalAnalysisforSemiCompRisks")
```

## Usage
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
simulation(2, 'null', sample_size = 1000, repeat_size = 100, get_variance = 'asymptotic')
simulation(2, 'alter', sample_size = 1000, repeat_size = 100, get_variance = 'asymptotic')
```
### REVEAL-HBV
```r
library(CausalAnalysisforSemiCompRisks)
# The result of analyzing the REVEAL-HBV dataset as discussed in Section 7 in our paper.
# We set downsample as 1 and get_variance as c('aymptotic', 'bootstrap') in our paper.
result = CASCR(REVEAL_HBV, downsample = 5, get_variance = 'asymptotic', plot_result = T)

# To get a quick result including bootstrap variance, one can run the following code.
# It will take around 15 minutes.
result = CASCR(REVEAL_HBV, downsample = 60, get_variance = c('asymptotic', 'bootstrap'), plot_result = T)

# The result of sensitivity analysis as presented in Supplement Material in Section 11.
result = CASCR(REVEAL_HBV, sen_ana = T, get_variance = NULL)
```

## More details
More details such as the input data type, parallel programming, and the output can be found using
```r
?CASCR
```
