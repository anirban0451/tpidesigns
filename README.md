
tpidesigns Package
================
Anirban Chakraborty

  - [tpidesigns: Dose Escalation Decisions based on Topxicity Probability Interval](#tpidesigns-introduction)
  - [Installation](#installation)
  - [Usage](#usage)
      - [1. Updation of Posterior Paramaters](#Updation-of-Posterior-Paramaters)
      - [2. Calculation of Unit Probability Mass](#Calculation-of-Unit-Probability-Mass)
      - [3. Decision Making Based on the number of DLT's](#Decision-Making-Based-on-the-number-of-DLT's)
      - [4. Graphical Plot of Unit Probability Mass](#Graphical-Plot-of-Unit-Probability-Mass)
      - [5. A Tabular Display of Dose Escalation Decisions](#A-Tabular-Display-of-Dose-Escalation-Decisions)
      

## tpidesigns: Introduction

This package is developed to build decisions based on Toxicity Probability Interval Designs which are very much familiar in Oncology early phase Clinical trials. This package primarily focuses on three designs - TPI, mTPI and mTPI-2. These Decisions are mostly driven by the amount of toxicity proportion suggested by the Clinician for a sample, the Prior Beta distribution (simple/ mixture) and several thresholds. Based on the importance of the particular drug, the Prior Distribution of the Toxicity distribution can change to simple or mixture.For documentation and conformability purpose, Toxicity will be termed as Dose Limiting Toxicity or DLT as per Oncology standard (Definition - [Here](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/dose-limiting))

By building this package, the user is provided with flexibilities to handle the parameters accordingly to find proper decisions in different scenarios. 

## Installation

``` r
devtools::install_github("anirban0451/tpidesigns")
```

## Usage

``` r
library(tpidesigns)
```

### 1\. Updation of Posterior Paramaters

In general, when a person (patient) is treated in Clinical Trials, his likelihood of developing DLT follows a Bernoulli Distribution and for a sample of more than one people, this likelihood will follow Binomial Distribution (Each patients independently may or may not develop DLT, so the number of DLT ' s in a sample will follow Binomial distribution, [Learn](https://en.wikipedia.org/wiki/Bernoulli_distribution#Related_distributions) more). `weights_formulate` takes into account the Binomial likelihood and the Prior distribution and hence develop the posterior distribution paramaters. Essentially, it generates a list of those parameters.

``` r
wt = runif(1) #generates a random weight on the first part of the prior distribution
params_updated = weights_formulate(w = wt, x = 1, n = 3, a1 = 1, b1= 1, a2 = 3, b2 = 5)
params_updated
```
Essentially,`w, a1, a2, b1, b2` are very crucial parameters for inputting the priors. In special scenarios, when there is only one prior, w should be input as either of 0 or 1, and only a1 and b1 should be input. Otherwise, a warning will show up.

### 2\. Calculation of Unit Probability Mass

The function `UPM` calculates Unit Probability Mass (Refer to [Ji et al(2010)](https://journals.sagepub.com/doi/pdf/10.1177/1740774510382799)) for a particular Interval within the range of Beta distribution.

``` r
#generating two random points within (0,1) for the Interval (not important, just for the sake of example)
interval = runif(2, min = 0 , max = 1)
a_min =  min(interval); b_max = max(interval)
#Generating a weight value
wt = runif(1) 
UPM(w = wt, a = a_min, b = b_max, a1 = 3, b1 = 6, a2 = 2, b2 = 5)
```
It is very important to note that a must be the minimum value in the interval and b must be the maximum value, otherwise the function will not work.
### 3\. Decision Making Based on the number of DLT's

When we have information about the number of DLT 's from a sample and we also have information about the prior distribution, then based on the target toxicity proportion (This a proportion of number of DLT 's in the sample assumed by mostly the Clinicians, he says the Dose is working properly when the Posterior mean of DLT proportion equates the assumed proportion) `decisiontpi` takes into account the decision making algorithm for Toxicity Probability Interval designs and gives us the decision for the particular sample. The decisions given are `Escalate (E)`, `Stay(S)`, `De-escalate(D)` and `Unacceptable(DU)`. Sample size must be at least 3.

``` r

n = 10 # sample size
pt = runif(1, min = 0.25, max = 0.35) #Target toxicity proportion 
x = sample.int(n, 1) #generating number of DLT 's in a sample
#Generating a weight value
wt = runif(1) 
decisiontpi(x = x, n = n, design = "mmtpi", pt = pt, e1 = 0.06, e2 = 0.04, eta = 0.95,
w = wt, a1 = 4, b1 = 3, a2 = 1, b2 = 1)
```
### 4\. Graphical Plot of DLT proportion

When it is confirmed that the Dose in not unacceptable, i.e. the function `decisiontpi` does not return the value `DU`, then one may want to understand the Posterior distribution of the Toxicity Proportion. `upmplot` takes into account the Decision Theoretic framework and develops the plot of Posterior Distribution of Toxicity Proportion (or DLT proportion). It also plots the values of Unit Probability Mass in different Toxicity Intervals (explained in the documentation). Sample size must be at least 3

```r
#Simulating some paramaters needed for the plot
n = 13 #must be a value >= 3
x = sample.int(n, 1)
pt = runif(1, min = 0.25, max = 0.35)
wt = runif(1)

require(ggplot2) #will be imported along with the main package
#Plotting of Posterior Distribution and UPM for mTPI-2(encoded as "mmtpi" design) design
upmplot(x = x, n = n, pt = pt, design = "mmtpi", w = wt, a1 = 1, a2 = 1, b1 = 4, b2 = 6) 

#Plotting of Posterior Distribution and UPM for mTPI design
upmplot(x = x, n = n, pt = pt, design = "mtpi", w = wt, a1 = 1, a2 = 1, b1 = 4, b2 = 6) 

```
### 5\. A Tabular Display of Dose Escalation Decisions

A person may want to know behaviour of a Drug with the help help of Prior information (encoded in the form of Prior Distribution in the parameters) and sample size, for example, how many DLT 's will allow the Clinician to increase the level of Drug, or how many DLT 's will conclude the current Drug level as unacceptably toxic. `tpitable` gives output for this scenario. When the maximum sample size (>= 3) and Prior information (w, a1, b1, a2, b2) is passed in the code, `tpitable` gives us a table where number of DLT's corresponding to the Decisions (explained in #3.) (more explanation of the Thresholds may be found in the documentation)

```r
#Simulating some paramaters needed for the table
nmax = 13 #must be a value >= 3
pt = runif(1, min = 0.25, max = 0.35)
wt = runif(1)

#Table for mTPI-2(encoded as "mmtpi" design) design
tpitable(nmax = nmax, pt = pt, eta = 0.95, design = "mmtpi", w = wt, a1 = 1, a2 = 1, b1 = 4, b2 = 6) 

#Table for mTPI design
tpitable(nmax = nmax, pt = pt, eta = 0.95, design = "mtpi", w = wt, a1 = 1, a2 = 1, b1 = 4, b2 = 6) 
```
