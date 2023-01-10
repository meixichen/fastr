
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpatialGEV

<!-- badges: start -->
<!--
[![R-CMD-check](https://github.com/meixichen/fastr/workflows/R-CMD-check/badge.svg)](https://github.com/meixichen/fastr/actions)
-->
<!-- badges: end -->

*Meixi Chen, Martin Lysy*

------------------------------------------------------------------------

## Description

The ***fastr*** (**F**actor **A**nalysis of **S**pike **T**rains in
**R**) package provides methods for analyzing neural spike trains
simultaneously recorded from multiple neurons, with a focus on exploring
the correlation between neurons.

## Installation

Before installing ***fastr***, make sure you have ***TMB*** installed
following the instructions
[here](https://github.com/kaskr/adcomp/wiki/Download).

To download the development version of this package, run

``` r
devtools::install_github("meixichen/fastr")
```

## Factor analysis of multi-neuron spike trains

The main method in this package models the multi-neuron spike trains
within a factor analysis framework. Simply put, it assumes that the
unobserved neural dynamics
![\\boldsymbol{x}(t)\\in\\mathbb{R}^q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bx%7D%28t%29%5Cin%5Cmathbb%7BR%7D%5Eq "\boldsymbol{x}(t)\in\mathbb{R}^q"),
which determines the observed spike trains
![\\boldsymbol{y}(t)\\in\\mathbb{R}^q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7By%7D%28t%29%5Cin%5Cmathbb%7BR%7D%5Eq "\boldsymbol{y}(t)\in\mathbb{R}^q"),
can be written as
![\\boldsymbol{x}(t)=\\boldsymbol{\\Lambda}\\boldsymbol{f}(t) + \\boldsymbol{\\epsilon}(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bx%7D%28t%29%3D%5Cboldsymbol%7B%5CLambda%7D%5Cboldsymbol%7Bf%7D%28t%29%20%2B%20%5Cboldsymbol%7B%5Cepsilon%7D%28t%29 "\boldsymbol{x}(t)=\boldsymbol{\Lambda}\boldsymbol{f}(t) + \boldsymbol{\epsilon}(t)"),
where
![\\boldsymbol{f}(t)\\in\\mathbb{R}^d, \\ d\\ll q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bf%7D%28t%29%5Cin%5Cmathbb%7BR%7D%5Ed%2C%20%5C%20d%5Cll%20q "\boldsymbol{f}(t)\in\mathbb{R}^d, \ d\ll q")
contains the independent “factors” and
![\\boldsymbol{\\Lambda}\_{q\\times d}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5CLambda%7D_%7Bq%5Ctimes%20d%7D "\boldsymbol{\Lambda}_{q\times d}")
is the factor loading matrix. Model fitting is carried out using the
main functionality in this package: `fastr_fit()`.

## Usage

More to come here.
