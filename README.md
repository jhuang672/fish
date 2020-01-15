# fish example for stochastic simulators: overview with opportunities

This repository contains all code and documentation with Rmarkdown guides to the fish agent based model example in the stochastic simulation review paper. 

### Rmarkdown files in this directory: 

### Main results in the paper: 

* **fish_fits.md**: [surrogate model fits and ABC calibration results](https://github.com/jhuang672/fish/blob/master/fish_fits.md) in the paper, including homGP, hetGP, sequential hetGP surrogates, and ABC calibration, knitted from **fish_fits.Rmd**.
* **fish_fits.Rmd**: for homGP, hetGP, and sequential hetGP surrogates, and ABC calibration code. 

### Data generation process: 

* **fish_sim.Rmd**: for grided on-shot space-filling design simulation
* **fish_seq.Rmd**: for sequential design simulation using IMSPE criteria 

Note: data simulation processes using fish model requires to wrap up NetLogo code into R environment 
(work for both Linux and windows systems now, but there is a rJava issue with Mac). 
Both Rmarkdown files set up the environment first and wrap up NetLogo from R environment. Then, generated data is further wrapped into R and saved into csv files. 

### The **data** directory contains the simulated fish data: 

* **GridData.csv**: one-shot simulation data from grided space-filling desgin 
* **SeqData.csv**: sequential design data using the IMSPE criteria 

### Who do I talk to? ###

* The project was originated and is actively maintained by Jiangeng Huang <huangj@ucsc.edu>
