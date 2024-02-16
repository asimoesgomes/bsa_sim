# Cost-benefit analysis of an AMC for broad spectrum antivirals

## Acknowledgements
This repo contains the scripts used in the cost benefit analysis of broad spectrum antivirals developed by SPRIND for the MSA Innovation Challenge
It builds on https://github.com/wwiecek/covstretch

## Notes

In this analysis, we project the potential reduction in deaths and infections from the distribution of broad spectrum antivirals to infected individuals in a COVID-like pandemic.

We also compare various starting points for the distribution in relation to the beginning of the epidemic (immediate, 100 days, 1 year), and different efficacies against transmission and death.

All strategies are compared according to the total number of infections and deaths.

## Setup

All results can be replicated through the [run-all-cases.R](run-all-cases.R) script. This should generate all of the figures and underlying data structures that are presented in the deliverable 1.

The main scripts called by [run-all-cases.R](run-all-cases.R) are: 

* **[project_setup.R](project_setup.R)**: Loads epidemiological models (ordinary differential equations implemented with odin), auxiliary functions, and initialize parameters.
* **[prep-results.R](cases/prep-results.R)**: Defines function to run the simulations and estimates the model using the status quo strategy (only vaccines and no antivirals)
* **[general-example.R](cases/general-example.R)**: Generates results and figures for outcomes under the status quo strategy.
* **[general-example-and-main.R](cases/general-example-and-main.R)**: Generates death and infection curves comparing the status quo and the three scenarios with antivirals (varying distribution timing).
* **[lower_efficacy_baseline_grid.R](cases/lower_efficacy_baseline_grid.R)**: Generates comparison matrix of the status quo and the scenarios with antivirals (varying both timing and efficacy assumptions).

