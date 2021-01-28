# Fish example for stochastic simulators: overview with opportunities

This repository contains all code and documentation with Rmarkdown guides to the fish agent based model example in the stochastic simulation review paper. 

### What is in this directory? 

### Main results in the paper: 

* **fish_fits.md**: [surrogate model fits and ABC calibration results](https://github.com/jhuang672/fish/blob/master/fish_fits.md) in the paper, including homGP, hetGP, QK, sequential hetGP surrogates, and ABC calibration, knitted from **fish_fits.Rmd**. The model fits are in square root scale of the original count data. Both plots in square root and original scales are provided. In additon to these limited sized design and model fits, dense grided space-filling design plots reflecting the underlying dynamic "truth" are also provided.
* **fish_fits.Rmd**: code for homGP, hetGP, QK, and sequential hetGP surrogates, and ABC calibration for fish example. The models are fitted in square root of count data and predictions are made by transforming back to the original scale. 

### Data generation process: 

R code for data genenration: 
* **fish_sim.Rmd**: for grided on-shot space-filling design simulation.
* **fish_seq.Rmd**: for sequential design simulation using IMSPE criteria. 
* **fish_sim_2.Rmd**: for dense grided one-shot space-filling design simulation.

NetLogo program and R wrapper: 
* **code/Fish.nlogo**: the NetLogo program for fish simulation in NetLogo version 6.0.4
* **code/netlogo_fish_functions.R**: R wrapper for NetLogo

Note: data simulation processes using fish model requires wrapping up NetLogo code into an R environment.
[NetLogo 6.0.4](https://ccl.northwestern.edu/netlogo/index.shtml) version was used for this paper, which works for both Linux and Windows systems, but with an rJava issue in macOS system. 
Both Rmarkdown files set up the parameters first and call NetLogo to run from R environment. Then, generated data is further wrapped into R and saved into csv files. 

### The **data** directory contains the simulated fish data: 

* **data/GridData.csv**: one-shot simulation data from grided space-filling desgin. 
* **data/SeqData.csv**: sequential design data using the IMSPE criteria. 
* **data/GridData_2.csv**: one-shot simulation data from a dense grided space-filling design, the "truth" simulation.  

### Who do I talk to? ###

* The project was originated and is actively maintained by Jiangeng Huang <huangj@vt.edu>
