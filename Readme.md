<div align="center">

# SEITR Network Analysis

[![R](https://img.shields.io/badge/R-≥4.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Academic_Use-lightgrey)]()
[![Status](https://img.shields.io/badge/Status-Pilot_Complete-brightgreen)]()
[![Topologies](https://img.shields.io/badge/Networks-ER_·_BA_·_WS-orange)]()

</div>

---


This project implements a SEITR (Susceptible, Exposed, Infected, Treated, Recovered) model for network analysis using R. The model simulates the spread of a disease through a network of individuals, allowing for various network types and parameters.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [SEITR_network](#seitr_network)
    - [Parameters](#parameters)
  - [compare_experiment_sets](#compare_experiment_sets)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

---

## Installation

To use this project, you need to have R installed on your system. You can install the required packages using the following commands:

```r
install.packages(c("igraph", "deSolve", "ggplot2", "dplyr"))
```

### Installing the Package from GitHub
You can install the SEITR Network Analysis package directly from GitHub using the devtools package. First, make sure you have devtools installed:
```r
install.packages("devtools")
```
Then, use the ```install_github``` function to install the package:
```r
devtools::install_github("skaraoglu/SEITRNet")
```
### Loading the Package
After installing the package, you can load it using the library function:
```r
library(SEITRNet)
```
---

## Usage

### SEITR_network
The SEITR_network function performs SEITR network analysis.

#### Parameters

- n: (int) Number of nodes in the network.
- network_type: (str) The type of network to create. This can be one of the following:

  - ER: Erdős-Renyi random graph
  - BA: Barabasi-Albert scale-free network
  - WS: Watts-Strogatz small-world network
  - LN: Lattice network
  - RR: Random regular network

- n_par1: (float) First network parameter. Assigned for "p" argument in Erdős-Renyi and Watts-Strogatz graphs, assigned for "k" argument for Lattice and Random regular, n * n_par1 is used as "m" for Barabási-Albert networks.

- n_par2: (float) Second network parameter that is only needed for Watts-Strogatz networks to assign "k" argument.

- Lambda: (float) Birth rate.

- alpha1: (float) Treatment rate
- alpha2: (float) Recovery rate from treatment
- delta_I: (float) Death rate from infection
- delta_T: (float) Death rate from treatment
- mu: (float) Natural death rate
- beta1: (float) The infection rate parameter.
- beta2: (float) The exposure rate parameter.
- beta3: (float) The recovery rate parameter.
- initial_statuses: (int) The initial status of the nodes in the network. Where each status can be one of the following:

    S: Susceptible, E: Exposed, I: Infected, T: Treatment, R: Recovered

- t: (int) Time period.
- num_exp: (int) The number of experiments to run.
- verbose: (bool) Verbose output.

```r
SEITR_network <- function(network_type="ER", n=100, n_par1=.9, n_par2=10, Lambda=1.1, beta1=.8, beta2=.18, beta3=.02, alpha1=.1, alpha2=.055, delta_I=.03, delta_T=.03, mu=.01, S=85, E=5, I=10, Tt=0, R=0, N=100, t=100, num_exp = 10, verbose = F, state = NULL, parameters = NULL) {
  # Function implementation
}
```
### Compare Experiment Sets
The compare_experiment_sets function compares the results of multiple experiment sets. Currently inactive.

---

## Examples
Here are some examples of how to use the functions in this project:

### Example 1: SEITR Network Analysis
```r
# Perform SEITR network analysis
ws_p.1_k20 <- SEITR_network("WS", n_par1=0.1, n_par2=20, num_exp = 3)
er_p.2 <- SEITR_network("ER", n_par1=0.2, num_exp = 3)
ba_m.75 <- SEITR_network("BA", n_par1=0.75, num_exp = 3)
```

### Example 2: Compare Experiment Sets
```r
# Compare the results of multiple experiment sets
compare_experiment_sets(list(ws_p.1_k20, er_p.2, ba_m.75))
```

---

## Citation

If you use this pipeline or build on this work, please cite: Karaoglu, S., Imran, M. & McKinney, B.A. Network-based SEITR epidemiological model with contact heterogeneity: comparison with homogeneous models for random, scale-free and small-world networks. Eur. Phys. J. Plus 140, 551 (2025). https://doi.org/10.1140/epjp/s13360-025-06481-z
