# Readme

This repository contains the data and code used for the analysis of the research project

__Inferring the proportion of undetected cholera infections from serological and clinical surveillance in an immunologically naive population__

by Flavio Finger, Joseph Lemaitre, Stanley Juin, Brendan Jackson, Sebastian Funk, Justin Lessler, Eric Mintz, Patrick Dely, Jacques Boncy and Andrew S Azman

which is currently under review. A preprint is available on MedrXiv: https://doi.org/10.1101/2023.11.01.23297461

This code and data repository is archived on Zenodo under doi [10.5281/zenodo.10063169](https://zenodo.org/doi/10.5281/zenodo.10063169).

## Citations

If you use this work or want to reference it, please cite our article as:

> _Inferring the proportion of undetected cholera infections from serological and clinical surveillance in an immunologically naive population_
> Flavio Finger, Joseph Lemaitre, Stanley Juin, Brendan Jackson, Sebastian Funk, Justin Lessler, Eric Mintz, Patrick Dely, Jacques Boncy, Andrew Azman
> medRxiv 2023.11.01.23297461;
> doi: https://doi.org/10.1101/2023.11.01.23297461


## Data

The incidence and serology datasets are located in the `data` folder:

The `incidence_clean.csv` file contains the clinical incidence data. Columns are:

- `date` : date of reporting
- `cas_vus` : reported cholera cases
- `cas_vus_lt5` : reported cholera cases < 5 years old
- `cas_vus_geq5` : reported cholera cases >= 5 years old

The `serology_clean.csv` file contains the serological data. Columns are:

- `date` : date serosamples taken
- `age` : age group (`2-4`, `>=5` or `NA` if missing)
- `vibriocidalMAX_titer` : measured vibriocidal titer value (maximum of Ogawa and Inaba) or `NA` if missing
- `n` : number of people sampled on `date` in age group `age` with titer value `vibriocidalMAX_titer`

See Jackson et al. (2013) [^1] for a detailed description of the data and survey methodology.

## Analysis code

The code is located in the `src` folder.

### Vibriocidal decay model

The vibriocidal decay model is implemented in Python, using PyMC3. The code is in the `src/py` folder

#### Requirements & versions

We are running our analysis in PyMC v4.4.0 and ArviZ v0.16.1. To run the code on any machine, the easiest way is to create an Anaconda environment as:

```bash
conda create -c conda-forge -n haiti-sero_pymc4 pymc=4 seaborn scipy arviz openpyxl jax numpyro pyreadr numpy pandas arviz click ipykernel pyreadr python=3.11
conda activate haiti-sero_pymc4
# Create a jupyter kernel with this environment
python -m ipykernel install --user --name haiti-sero_pymc4_env --display-name "Python (haiti-sero_pymc4)"
```

To reproduce the code using the exact package versions we used for the manuscript, we provide two conda environment files:
- `environment_exact.yml` has the exact package versions used for this project to fit and analyze the results. It may only work on Mac OS.
- `environment_crossplatform.yml` has just the packages that are required for installation, and is equivalent to the command line above. It should install on any platforms.

See [Anaconda's documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#sharing-an-environment) on conda environments.


#### Running the code

##### File descriptions

- `src/py/Fit-hist_AnalysisAndSetup.ipynb` holds the development version of the PyMC model, and the analysis code that produces all the figures in the main text and in the supplementary information.
- `src/py/calibration_fit-hist.py` is the calibration script, with the same model as the above notebook. It requires one system argument: the identifier of the model specification: `0` for aged 2-4, `1` for age 5+, and `2` for the full model with everyone.
- `src/py/batch_fit-hist.run` is a slurm file to run the above calibration script on an HPC cluster, with all model specifications at the same time.
- `src/py/utils.py` has the code to load, clean and prepare the two datafiles along with some helper functions.
- `src/py/priors from Jones et al/` contains the code that has been used to extract the priors from the data from Jones et al. (2022) [^2], these are used as the Multivariate peak titer priors.


The sensitivity analyses shown in the supplementary information are in the files:

- `src/py/calibration_fit-hist_rate_sensitivity.py` which is similar to `src/py/calibration_fit-hist.py`, except that it only models age 2-4 and the command line argument instead is a step of perturbation for the decay rate.
- `src/py/batch_fit-hist-sensitivity.run` runs the above script on an HPC cluster for different decay rates

The steps to reproduce our results involves:

1. Calibrating the model using `python src/py/calibration_fit-hist.py` with arguments 0,1 and 2 or doing all three at the same time with the batch runfile.
2. Run notebook `src/py/Fit-hist_AnalysisAndSetup.ipynb` which produces all figures and numerical output shown in the manuscript with their credible interval and the diagnostic checks of the model.

### Attack rate estimates

The code for the additional attack rate estimates is written in R and Stan.

#### Requirements & versions

The analyses were performed using R version 4.3.1.
For the Gaussian mixture model we used rstan version 2.32.3 with Stan version 2.26.1.

Packages required can be installed with
`install.packages(c("rmarkdown", "knitr", "here", "dplyr", "tidyr", "ggplot2", "magrittr", "rstan", "readr"))`

#### Running the code

Compare different attack rate and infection rate estimates:
`rmarkdown::render("src/Rmd/attack_rate_estimates.rmd")`

The output html file will be saved in the same folder.


## References

[^1]: _B. R. Jackson et al._, “Seroepidemiologic Survey of Epidemic Cholera in Haiti to Assess Spectrum of Illness and Risk Factors for Severe Disease,” The American Journal of Tropical Medicine and Hygiene, vol. 89, no. 4, pp. 654–664, Oct. 2013, doi: [10.4269/ajtmh.13-0208](https://doi.org/10.4269/ajtmh.13-0208).


[^2]: _F. K. Jones et al._, “Identifying Recent Cholera Infections Using a Multiplex Bead Serological Assay,” mBio, vol. 13, no. 6, pp. e01900-22, Dec. 2022, doi: [10.1128/mbio.01900-22](https://doi.org/10.1128/mbio.01900-22).