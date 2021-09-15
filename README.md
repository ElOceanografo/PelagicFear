# Analysis code for "Fear structures the ocean's pelagic zone"

[![DOI](https://zenodo.org/badge/397781245.svg)](https://zenodo.org/badge/latestdoi/397781245)

This directory contains the code necessary to reproduce the figures and analysis in Urmy, S.S. and K.J. Benoit-Bird (2021), "Fear structures the ocean's pelagic zone," https://doi.org/10.1016/j.cub.2021.09.003

## Data
The data required for these analyses are archived at Dryad: https://doi.org/10.5061/dryad.44j0zpcf1.  They will need to be downloaded before the code can be run.

## Required software

Analysis was performed using the open-source languages Julia (v1.6.2) and R (v4.0.5). These can be downloaded from http://julialang.org and http://r-project.org.

### Julia packages

All required Julia packages are recorded in the files `Project.toml` and `Manifest.toml`.  They can be installed automatically, at the exact versions used for the analysis, by opening a Julia REPL (read-evaluate-print loop, i.e. the Julia command line) in this directory and performing the following steps:

1. Type `]` to enter package manager mode. The prompt will change from `julia>` to `(@1.6) pkg>`.
2. Type `activate .` and hit Enter to activate the project's package environment. The prompt will change to `(Fear) pkg>`.
3. Enter `instantiate`. Julia will install and precompile all the required packages (may take a few minutes).  When done, you can press backspace to exit the package manager.

### R packages

The R scripts require the following packages.  The versions used are noted in parentheses:

* `dplyr` (v1.0.5)
* `tidyr` (v.1.3)
* `arrow` (v4.0.0)
* `lubridate` (v1.7.10)
* `ggplot2` (v3.3.3)
* `viridis` (v0.6.0)
* `scales` (v1.1.1)


## Running the analyses

The Julia script `school_detection.jl` should be run first.  If Julia is on your system path, you can run `julia school_detection.jl` on your system's command line. Alternatively, after doing steps 1-3 in the "Julia packages" section above, you can run `include("src/school_detection.jl")` from the Julia prompt, or run the script interactively through VSCode (https://www.julia-vscode.org/) or Juno (https://junolab.org/).

The R scripts can be run in any order, either from the system command line (e.g., `R CMD BATCH echogram_plots.R`), or interactively from the R GUI or RStudio.
