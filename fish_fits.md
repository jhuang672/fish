Fish Example in Stochastic Simulators
================
Jiangeng Huang
June 27, 2020

This file generates three types of surrogate models on the fish
agent-based model:

  - **HomGP**: homoscedastic Gaussian process surrogate on one-shot
    grided space-filling design
  - **HetGP**: heteroskedastic Gaussian process surrogate on one-shot
    grided space-filling design
  - **Sequential design**: heteroskedastic Gaussian process surrogate on
    sequential design based on IMSPE criteria

and calibration of the fish model using:

  - **ABC**: approximate Bayesian computation

A square root transformation has been taken to all the responses (which
are count data) in the model fitting step, in order to better fit the
Gaussian error assumptions. Both plots in square root and original
scales are shown.

For grided design and sequential design, both data generation processes
require to wrap up NetLogo through R environment, see the following two
Rmarkdown files in this directory on data geneeration:

  - **fish\_sim.Rmd** for grided on-shot space-filling design simulation
  - **fish\_seq.Rmd** for sequential design simulation using IMSPE
    criteria

-----

-----

Fitting standard homoscedastic and heteroscedastic GP surrogates using
the hetGP package:

``` r
library(hetGP)
```

## HomGP

We can fit the homGP emulator, and then make predictions

Step 1: read in the CSV file of fish simulator runs and prepare the data

``` r
fish <- read.csv("data/GridData.csv")

## set up boundary for X values: 
lower_b <- 150
upper_b <- 4000

## Scale X in [150, 4000] into [0 ,1] unit cube
X <- ( fish$population_size - lower_b )/ (upper_b - lower_b)  
Y <- sqrt(fish$recaptures)
```

Step 2: Train the model

``` r
mod.hom <- mleHomGP(X = X, Z = Y, lower = 0.0001, upper = 10)
```

Step 3: Make predictions

``` r
xgrid <- seq(0, 1, length = 1000) #values to make predictions
p.hom <- predict(mod.hom, matrix(xgrid, ncol = 1)) #make predictions
pvar.hom <- p.hom$sd2 + p.hom$nugs #get total preditive variance
```

Step 4: make plots

``` r
## in sqrt scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, p.hom$mean, type = "l",xlab = "Population", ylab = "Square Root of Number of Marked in Recapture",   
     main ="Homoskedastic Gaussian process surrogate on gridded design ", ylim = c(0, 10))
points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.05, p.hom$mean, sqrt(pvar.hom)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.95, p.hom$mean, sqrt(pvar.hom)), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20a-1.png)<!-- -->

``` r
## transform back in original scale for population: 
plot( xgrid* (upper_b - lower_b) + lower_b, (p.hom$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Homoskedastic Gaussian process surrogate on gridded design", ylim = c(0, 85))
points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.05, p.hom$mean, sqrt(pvar.hom))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.95, p.hom$mean, sqrt(pvar.hom))^2), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20a-2.png)<!-- -->

## HetGP

We can then repeat the whole process, but for hetGP. The only difference
here is in step 2.

Step 1: read in the CSV file of fish simulator runs and prepare the data

``` r
fish <- read.csv("data/GridData.csv")

## set up boundary for X values: 
lower_b <- 150
upper_b <- 4000

## Scale X in [150, 4000] into [0 ,1] unit cube
X <- ( fish$population_size - lower_b )/ (upper_b - lower_b)  
Y <- sqrt(fish$recaptures)
```

Step 2: Train the
model

``` r
mod.het <- mleHetGP(X = X, Z = Y, lower = 0.0001, upper = 10) #hetGP now instead
```

Step 3: Make predictions

``` r
xgrid <- seq(0, 1, length = 1000) #values to make predictions
p.het <- predict(mod.het, matrix(xgrid, ncol = 1)) #make predictions
pvar.het <- p.het$sd2 + p.het$nugs #get total preditive variance
```

Step 4: make plots

``` r
## in sqrt scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, p.het$mean, type = "l",xlab = "Population", ylab = "Square Root of Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on gridded design", ylim = c(0, 10))
points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.05, p.het$mean, sqrt(pvar.het)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.95, p.het$mean, sqrt(pvar.het)), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20b-1.png)<!-- -->

``` r
## transform back in original scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, (p.het$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on gridded design", ylim = c(0, 85))
points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.05, p.het$mean, sqrt(pvar.het))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.95, p.het$mean, sqrt(pvar.het))^2), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20b-2.png)<!-- -->

## Sequential Design

And again, but using the sequential design (and using hetGP). Not shown
here is the mechanism to implement the sequential design; this is
invariably linked to the simulator itself, but it is also fairly
straightforward. See the documentation for hetGP, or the Ocean Rmarkdown
file for more information. The only difference here is in Step 1.

Step 1: read in the CSV file of fish simulator runs and prepare the data

``` r
fish <- read.csv("data/SeqData.csv") #sequential data now instead

## set up boundary for X values: 
lower_b <- 150
upper_b <- 4000

## Scale X in [150, 4000] into [0 ,1] unit cube
X <- ( fish$population_size - lower_b )/ (upper_b - lower_b)  
Y <- sqrt(fish$recaptures)
```

Step 2: Train the model

``` r
mod.seq <- mleHetGP(X = X, Z = Y, lower = 0.0001, upper = 10)
```

Step 3: Make predictions

``` r
xgrid <- seq(0, 1, length = 1000) #values to make predictions
p.seq <- predict(mod.seq, matrix(xgrid, ncol = 1)) #make predictions
pvar.seq <- p.seq$sd2 + p.seq$nugs #get total preditive variance
```

Step 4: make plot

``` r
## in sqrt scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, p.seq$mean, type = "l",xlab = "Population", ylab = "Square Root of Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on sequential design", ylim = c(0, 10))
points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.05, p.seq$mean, sqrt(pvar.seq)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.95, p.seq$mean, sqrt(pvar.seq)), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20c-1.png)<!-- -->

``` r
## transform back in original scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, (p.seq$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on sequential design", ylim = c(0, 85))
points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.05, p.seq$mean, sqrt(pvar.seq))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.95, p.seq$mean, sqrt(pvar.seq))^2), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20c-2.png)<!-- -->

## ABC

Applying ABC to one of these surrogates is very easy.

Step 1: we \`\`import’’ our observed data

``` r
# observed data
yobs <- sqrt(25) 
```

Step 2: we simulate from our surrogate many times. Here we choose to use
the standard hetGP surrogate, although any could be used.

``` r
set.seed(123)
N <- 100000 ## number of ABC draws
xsim <- seq(0, 1, length.out = N) ## select x-values to try
p.het <- predict(mod.het, matrix(xsim, ncol = 1)) ## predictions using hetGP
pvar.het <- p.het$sd2 + p.het$nugs ## prediction var

ysim <- rnorm(N, mean = p.het$mean, sd = sqrt(pvar.het)) ## simulate from normal (hetGP gives mean and variance of this normal)
ysim <- as.integer(ysim) ## take the integer part (can't have half a fish...)
```

Step 3: we check which simulated values equal the observed value

``` r
## ABC
x.ABC <- xsim[which(ysim == yobs)] ## check which simulations equal the observation
x.ABC.native <- x.ABC * (upper_b - lower_b) + 150 ## convert these matching x values into the right scale (not (0,1))
```

Step 4: plot the resulting posterior histogram

``` r
## plot
hist(x.ABC.native, breaks = 30, probability = T, main="ABC posterior using hetGP surrogate", xlim = c(0, 1500), xlab = 'Population') ## plot histogram
```

![](fish_fits_files/figure-gfm/ABC%20plotting-1.png)<!-- -->

``` r
#lines(density(x.ABC.native), col = 'red', lwd = 2) ##add densitiy estimate (optional)
```

Doing this without a surrogate, which is not recommended, would simply
involve replacing step 2 with simulating from the simulator (rather than
from the surrogate)
