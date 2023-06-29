
<!-- badges: start -->

[![R-CMD-check](https://github.com/soylentOrange/sdarr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/soylentOrange/sdarr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdarr

The sdarr-package provides a R implementation of the SDAR algorithm for
Slope Determination by Analysis of Residuals as standardized in [ASTM
E3076-18](https://doi.org/10.1520/E3076-18).  
It allows for automated and objective linear-fitting of mechanical
test-data. See a detailed description of the original algorithm in the
[NIST Technical Note 2050 by E.
Lucon](https://doi.org/10.6028/NIST.TN.2050).

As the original algorithm heavily uses linear regressions, a faster
version was implemented here, which finds the optimum region for linear
fitting by random sub-sampling.

## Installation

You can install the development version of sdarr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("soylentOrange/sdarr")
```

## Acknowledgements

The sdarr-package was created for the analysis of mechanical test data
for the project LOBio. The project LOBio is funded by the German
ministry of education and research (BMBF) under grant number 13XP5174C.

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(sdarr) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
