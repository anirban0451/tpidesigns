
tpidesigns Package
================
Anirban Chakraborty

  - [tpidesigns: Dose Escalation Decisions based on Topxicity Probability Interval](#tpidesigns-introduction)
  - [Installation](#installation)
  - [Usage](#usage)
      - [1. Updation of Posterior Paramaters](#Updation-of-Posterior-Paramaters)
      - [2. Calculation of Unit Probability Mass](#Calculation-of-Unit-Probability-Mass)
      - [3. Decision Making Based on the number of DLT's](#Decision-Making-Based-on-the-number-of-DLT's)
      

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

This function calculates Unit Probability Mass (Refer to [Ji et al(2010)](https://journals.sagepub.com/doi/pdf/10.1177/1740774510382799)) for a particular Interval within the range of Beta distribution.

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

