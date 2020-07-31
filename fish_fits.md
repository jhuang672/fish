Fish Examples in Stochastic Simulators
================
Jiangeng Huang
July 15, 2020

This file generates three types of surrogate models on the fish
agent-based model:

  - **HomGP**: homoscedastic Gaussian process surrogate on one-shot
    grided space-filling design.
  - **HetGP**: heteroskedastic Gaussian process surrogate on one-shot
    grided space-filling design.
  - **QK**: quantile kriging surrogate on one-shot grided space-filling
    design.
  - **Sequential design**: heteroskedastic Gaussian process surrogate on
    sequential design based on IMSPE criteria.
  - **Superimposed plots**: superimposed plots of one-shot homGP,
    one-shot hetGP, sequential hetGP fits over the “truth” quantiles on
    dense grided design.

and calibration of the fish model using:

  - **ABC**: approximate Bayesian computation.

A square root transformation has been taken to all the responses (which
are count data) in the model fitting step, in order to better fit the
Gaussian error assumptions. Both plots in square root and original
scales are shown.

For grided design and sequential design, both data generation processes
require to wrap up NetLogo through R environment, see the following two
Rmarkdown files in this directory on data geneeration:

  - **fish\_sim.Rmd** for grided on-shot space-filling design
    simulation.
  - **fish\_seq.Rmd** for sequential design simulation using IMSPE
    criteria.
  - **fish\_sim\_2.Rmd** for dense grided one-shot space-filling design
    simulation.

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
     main ="Homoskedastic Gaussian process surrogate on gridded design ", ylim = c(0, 10), col = 2)
points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.025, p.hom$mean, sqrt(pvar.hom)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.975, p.hom$mean, sqrt(pvar.hom)), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20a-1.png)<!-- -->

``` r
## transform back in original scale for population: 
plot( xgrid* (upper_b - lower_b) + lower_b, (p.hom$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Homoskedastic Gaussian process surrogate on gridded design", ylim = c(0, 85), col = 2)
points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.025, p.hom$mean, sqrt(pvar.hom))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.975, p.hom$mean, sqrt(pvar.hom))^2), col = 2, lty = 2)
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

Step 2: Train the model

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
     main ="Heteroskedastic Gaussian process surrogate on gridded design", ylim = c(0, 10), col = 2)
points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.025, p.het$mean, sqrt(pvar.het)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.975, p.het$mean, sqrt(pvar.het)), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20b-1.png)<!-- -->

``` r
## transform back in original scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, (p.het$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on gridded design", ylim = c(0, 85), col = 2)
points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.025, p.het$mean, sqrt(pvar.het))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.975, p.het$mean, sqrt(pvar.het))^2), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20b-2.png)<!-- -->

## Quantile Kriging

We can also fit a quantile kriging emulator

Step 1: read in the CSV file of fish simulator runs and prepare the
data. For QK, this involves a fair bit of re-arranging.

``` r
fish <- read.csv("data/GridData.csv")

## set up boundary for X values: 
lower_b <- 150
upper_b <- 4000

## Scale X in [150, 4000] into [0 ,1] unit cube
X <- ( fish$population_size - lower_b )/ (upper_b - lower_b)  
Y <- fish$recaptures

#need to re-arrange the data
recaps_fit <- matrix(fish$recaptures, ncol = 20, nrow = 20, byrow = TRUE) #first arrange output by unique input
recaps_fit_qs <- apply(recaps_fit, 2, quantile, probs=c(0.05, 0.275, 0.5, 0.725, 0.95)) #then calculate empirical quantiles
recaps_fit_qs_y <- as.numeric(recaps_fit_qs) #then flatten down into a 1D vector
recaps_fit_qs_x <- rep(unique(X),each = 5) #and get the input values to include q as an input
recaps_fit_qs_q <- rep(c(0.05, 0.275, 0.5, 0.725, 0.95), 20)
```

Step 2: Train the model

``` r
mod.QK <- mleHomGP(X = cbind(recaps_fit_qs_x, recaps_fit_qs_q), Z = recaps_fit_qs_y)
```

Step 3: Make predictions

``` r
xgrid <- seq(0, 1, length = 1000) #values to make predictions

#and make predictions (need to predict for each quantile we are interested in)
p.QK.05 <- predict(mod.QK, cbind(xgrid, 0.05) )
p.QK.25 <- predict(mod.QK, cbind(xgrid, 0.25) )
p.QK.5 <- predict(mod.QK, cbind(xgrid, 0.5) )
p.QK.75 <- predict(mod.QK, cbind(xgrid, 0.75) )
p.QK.95 <- predict(mod.QK, cbind(xgrid, 0.95) )
```

Step 4: make plots

``` r
plot(xgrid* (upper_b - lower_b) + lower_b, p.QK.05$mean, type = 'l', lty = 2, col=c("red"), ylim = c(0, 85),  xlab = "Population", ylab = "Number of Marked in Recapture", main ="Quantile kriging surrogate on gridded design ")
points(fish$population_size , fish$recaptures)
lines(xgrid* (upper_b - lower_b) + lower_b, p.QK.25$mean, type = 'l', lty = 2, col=c("blue"))
lines(xgrid* (upper_b - lower_b) + lower_b, p.QK.5$mean, type = 'l', lty = 1, col=c("purple"))
lines(xgrid* (upper_b - lower_b) + lower_b, p.QK.75$mean, type = 'l', lty = 2, col=c("blue"))
lines(xgrid* (upper_b - lower_b) + lower_b, p.QK.95$mean, type = 'l', lty = 2, col=c("red"))
legend('topright', c(paste0(c(95, 75, 50, 25, 5), "%")), lty = c(2,2,1,2,2, 2), col=c("red", "blue", "purple", "blue", "red"))
```

![](fish_fits_files/figure-gfm/make%20plots%20QK-1.png)<!-- -->

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
     main ="Heteroskedastic Gaussian process surrogate on sequential design", ylim = c(0, 10), col = 2)
points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.025, p.seq$mean, sqrt(pvar.seq)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.975, p.seq$mean, sqrt(pvar.seq)), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20c-1.png)<!-- -->

``` r
## transform back in original scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, (p.seq$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on sequential design", ylim = c(0, 85), col = 2)
points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.025, p.seq$mean, sqrt(pvar.seq))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.975, p.seq$mean, sqrt(pvar.seq))^2), col = 2, lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20c-2.png)<!-- -->

## Obtaining the “Truth”:

500 dense replications at the exactly same 20 grided locations of the
one-shot space-filling design.

``` r
## Read in "truth" data:
fish_2 <- read.csv("data/GridData_2.csv")

## Re-organize for quantiles plots: 
recaps <- matrix(fish_2$recaptures, ncol = 20, nrow = 500, byrow = TRUE)
```

And plot it:

``` r
## Plots of the raw truth simulation: 
plot(fish_2$population_size, fish_2$recaptures, type = "p",xlab = "Population", ylab = "Number of Marked in Recapture", ylim = c(0, 85))
```

![](fish_fits_files/figure-gfm/make%20plots%20d-1.png)<!-- -->

``` r
boxplot(recaps, main = "Boxplots of the 'truth' simulations")   
```

![](fish_fits_files/figure-gfm/make%20plots%20d-2.png)<!-- -->

``` r
## Quantiles plots: 
recaps_qs <- apply(recaps, 2, quantile, probs=c(0.025, 0.5, 0.975))

## Quantile line plots
plot(fish_2$population_size[1:20], recaps_qs[2,], ylim = c(0, 85),  xlab = "Population", ylab = "Number of Marked in Recapture", type = "l", main = "Quantiles of 'truth' simulations")
lines(fish_2$population_size[1:20], recaps_qs[1,], lty = 2)
lines(fish_2$population_size[1:20], recaps_qs[3,], lty = 2)
```

![](fish_fits_files/figure-gfm/make%20plots%20d-3.png)<!-- -->

## Superimposed plots

### HomGP with “truth”:

``` r
## in sqrt scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, p.hom$mean, type = "l",xlab = "Population", ylab = "Square Root of Number of Marked in Recapture",   
     main ="Homoskedastic Gaussian process surrogate on gridded design", ylim = c(0, 10), col = 2)
#points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.025, p.hom$mean, sqrt(pvar.hom)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.975, p.hom$mean, sqrt(pvar.hom)), col = 2, lty = 2)

## plot pointwise quantiles only: 
lines(fish_2$population_size[1:20], sqrt(recaps_qs[2,]), col = 1)
lines(fish_2$population_size[1:20], sqrt(recaps_qs[1,]), col = 1, lty = 2)
lines(fish_2$population_size[1:20], sqrt(recaps_qs[3,]), col = 1, lty = 2)
legend("topright", c("homGP", "Truth"), pch = c(NA, NA, NA), lty = c(1, 1), col=c(2, 1))
```

![](fish_fits_files/figure-gfm/make%20plots%20e-1.png)<!-- -->

``` r
## In original units: 
## transform back in original scale for population: 
plot( xgrid* (upper_b - lower_b) + lower_b, (p.hom$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Homoskedastic Gaussian process surrogate on gridded design", ylim = c(0, 85), col = 2)
#points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.025, p.hom$mean, sqrt(pvar.hom))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.975, p.hom$mean, sqrt(pvar.hom))^2), col = 2, lty = 2)


## plot pointwise quantiles only: 
lines(fish_2$population_size[1:20], recaps_qs[2,], col = 1)
lines(fish_2$population_size[1:20], recaps_qs[1,], col = 1, lty = 2)
lines(fish_2$population_size[1:20], recaps_qs[3,], col = 1, lty = 2)
legend("topright", c("homGP", "Truth"), pch = c(NA, NA, NA), lty = c(1, 1), col=c(2, 1))
```

![](fish_fits_files/figure-gfm/make%20plots%20e-2.png)<!-- -->

### HetGP with “truth”:

``` r
## in sqrt scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, p.het$mean, type = "l",xlab = "Population", ylab = "Square Root of Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on gridded design", ylim = c(0, 10), col = 2)
#points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.025, p.het$mean, sqrt(pvar.het)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.975, p.het$mean, sqrt(pvar.het)), col = 2, lty = 2)


## plot pointwise quantiles only: 
lines(fish_2$population_size[1:20], sqrt(recaps_qs[2,]), col = 1)
lines(fish_2$population_size[1:20], sqrt(recaps_qs[1,]), col = 1, lty = 2)
lines(fish_2$population_size[1:20], sqrt(recaps_qs[3,]), col = 1, lty = 2)
legend("topright", c("hetGP", "Truth"), pch = c(NA, NA, NA), lty = c(1, 1), col=c(2, 1))
```

![](fish_fits_files/figure-gfm/make%20plots%20f-1.png)<!-- -->

``` r
## In original units: 
## transform back in original scale for population: 
plot( xgrid* (upper_b - lower_b) + lower_b, (p.het$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on gridded design", ylim = c(0, 85), col = 2)
#points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.025, p.het$mean, sqrt(pvar.het))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.975, p.het$mean, sqrt(pvar.het))^2), col = 2, lty = 2)


## plot pointwise quantiles only: 
lines(fish_2$population_size[1:20], recaps_qs[2,], col = 1)
lines(fish_2$population_size[1:20], recaps_qs[1,], col = 1, lty = 2)
lines(fish_2$population_size[1:20], recaps_qs[3,], col = 1, lty = 2)
legend("topright", c("hetGP", "Truth"), pch = c(NA, NA, NA), lty = c(1, 1), col=c(2, 1))
```

![](fish_fits_files/figure-gfm/make%20plots%20f-2.png)<!-- -->

### QK with “truth”:

``` r
recaps_qs_2 <- apply(recaps, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95)) #truth for different quantiles


plot(xgrid* (upper_b - lower_b) + lower_b, p.QK.05$mean, type = 'l', lty = 2, col=c("red"), ylim = c(0, 85),  xlab = "Population", ylab = "Number of Marked in Recapture", main ="Quantile kriging surrogate on gridded design ")
lines(xgrid* (upper_b - lower_b) + lower_b, p.QK.5$mean, type = 'l', lty = 1, col=c("purple"))
lines(xgrid* (upper_b - lower_b) + lower_b, p.QK.95$mean, type = 'l', lty = 2, col=c("red"))


## plot pointwise quantiles only: 
lines(fish_2$population_size[1:20], recaps_qs_2[3,], col = 1)
lines(fish_2$population_size[1:20], recaps_qs_2[1,], col = 1, lty = 2)
lines(fish_2$population_size[1:20], recaps_qs_2[5,], col = 1, lty = 2)
legend('topright', c(paste0(c(95, 50, 5), "%"), "Truth"), lty = c(2,1,2, 1), col=c("red", "purple", "red", "black"))
```

![](fish_fits_files/figure-gfm/make%20plots%20QK%202-1.png)<!-- -->

### Sequential hetGP with “truth”:

``` r
## in sqrt scale for population: 
plot(xgrid* (upper_b - lower_b) + lower_b, p.seq$mean, type = "l",xlab = "Population", ylab = "Square Root of Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on sequential design", ylim = c(0, 10), col = 2)
#points(fish$population_size , Y)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.025, p.seq$mean, sqrt(pvar.seq)), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, qnorm(0.975, p.seq$mean, sqrt(pvar.seq)), col = 2, lty = 2)

## plot pointwise quantiles only: 
lines(fish_2$population_size[1:20], sqrt(recaps_qs[2,]), col = 1)
lines(fish_2$population_size[1:20], sqrt(recaps_qs[1,]), col = 1, lty = 2)
lines(fish_2$population_size[1:20], sqrt(recaps_qs[3,]), col = 1, lty = 2)
legend("topright", c("SeqhetGP", "Truth"), pch = c(NA, NA, NA), lty = c(1, 1), col=c(2, 1))
```

![](fish_fits_files/figure-gfm/make%20plots%20g-1.png)<!-- -->

``` r
## In original units: 
## transform back in original scale for population: 
plot( xgrid* (upper_b - lower_b) + lower_b, (p.seq$mean)^2, type = "l",xlab = "Population", ylab = "Number of Marked in Recapture",   
     main ="Heteroskedastic Gaussian process surrogate on sequential design", ylim = c(0, 85), col = 2)
#points(fish$population_size , Y^2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.025, p.seq$mean, sqrt(pvar.seq))^2), col = 2, lty = 2)
lines(xgrid * (upper_b - lower_b) + lower_b, (qnorm(0.975, p.seq$mean, sqrt(pvar.seq))^2), col = 2, lty = 2)


## plot pointwise quantiles only: 
lines(fish_2$population_size[1:20], recaps_qs[2,], col = 1)
lines(fish_2$population_size[1:20], recaps_qs[1,], col = 1, lty = 2)
lines(fish_2$population_size[1:20], recaps_qs[3,], col = 1, lty = 2)
legend("topright", c("SeqhetGP", "Truth"), pch = c(NA, NA, NA), lty = c(1, 1), col=c(2, 1))
```

![](fish_fits_files/figure-gfm/make%20plots%20g-2.png)<!-- -->

## ABC

Applying ABC to one of these surrogates is very easy.

Step 1: we \`\`import’’ our observed data

``` r
# observed data
yobs <- 25 
```

Step 2: we simulate from our surrogate many times. Here we choose to use
the standard hetGP surrogate, although any could be used.

``` r
set.seed(123)
N <- 1000000 ## number of ABC draws
xsim <- seq(0, 1, length.out = N) ## select x-values to try
p.het <- predict(mod.het, matrix(xsim, ncol = 1)) ## predictions using hetGP
pvar.het <- p.het$sd2 + p.het$nugs ## prediction var

ysim <- rnorm(N, mean = p.het$mean, sd = sqrt(pvar.het)) ## simulate from normal (hetGP gives mean and variance of this normal)
ysim <- ysim^2 # undo sqrt transform
ysim <- as.integer(ysim) ## take the integer part (can't have half a fish...)
```

Step 3: we check which simulated values equal the observed value

``` r
## ABC
x.ABC <- xsim[which(ysim == yobs)] ## check which simulations equal the observation
x.ABC.native <- x.ABC * (upper_b - lower_b) + 150 ## convert these matching x values into the right scale (not (0,1))
print(length(x.ABC.native)) #how many accepted samples?
```

    ## [1] 3811

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
