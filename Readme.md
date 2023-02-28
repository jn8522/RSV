<h1 id="top"> Repo overview </h1>

-   [**src**](#section.src) --- source code – contains functions used by files in *scripts*

-   [**data**](#section.data) --- contains data used to fit the models as well as summary statistics on model fits

-   [**scripts**](#section.scripts) --- contains the executable code that does all the work

    -   **plots** --- contains code that generates the plots in
        **Repo/plots**

-   [**models**](#section.models) --- contains the fitted models

-   [**diagnostics**](#section.diagnostics) --- contains results from test runs of fitting the
    models

-   [**plots**](#section.plots) --- contains plots used in the paper

# Files

[<h2 id="section.src"> src </h2>](#top)

  -   *ODE.R* --- contains functions that implement the ODE models

  -   *fit.R* --- contains functions related to fitting the ODE models

  -   *plot.R* --- contains functions that help with plotting the model
        fits

[<h2 id="section.data"> data </h2>](#top)

-   **raw** --- contains raw data

    -   *bom_data.txt* --- contains daily weather observations for the
        Perth Airport weather station; this data cannot be shared under our agreements with the data providers

    -   *Flu_RSV.csv* --- list of hospital admissions and separations,
        including RSV test status; this data cannot be shared under our agreements with the data providers

    -   *school_terms.csv* --- list of term dates for Western Australian
        public schools

    -   *known.values.csv* --- list of values for model parameters that are
        fixed (i.e. are not being fit to the data)

-   **processed** --- contains data after being processed by
    **scripts**/*wrangling.R*

    -   *weather.school.csv* --- weekly observations of weather variables
        and whether school is in term; this data cannot be shared under our agreements with the data providers

    -   *RSV.csv* --- RSV case numbers by week; this data cannot be shared under our agreements with the data providers
- **statistics** --- contains summary statistics on model fits
	+ *cis.profile.base.RData* --- confidence intervals for the sinusoidal model calculated using profile likelihood
	+ *cis.profile.weather.Rdata* --- confidence intervals for the weather model calculated using profile likelihood
	+ *cis.wald.base.Rdata* --- Wald confidence intervals for the sinusoidal model
	+ *cis.wald.weather.Rdata* --- Wald confidence intervals for the weather model

[<h2 id="section.scripts"> scripts </h2>](#top)
These should be run in the order below.

-   *wrangling.R* --- processes the data from **data/raw** into the data in
    **data/processed**
- 	*fit.BackwardsSelection* --- fits different versions of the weather model using backwards selection to see which weather variables should be included in the final test against the base model. Models are compared using AIC and results stored as **models**/*fits.backwards.selection.Rdata*.
-   *fit.R* --- fits each model using 5-fold cross validation and
    stores the results as **models**/*fits.\_\_\_.CV.Rdata*; also fits each base model to the entire time series (this will help produce Figures 3-5) and stores the result as **models**/*fit.\_\_\_.Rdata*
-   *diagnostics.R* --- checks whether varying the initial values for
    the parameters affects the fit for either model; results stored in
    **diagnostics**/*diagnostic.fits.\_\_\_.Rdata*
-   *statistics.R* --- calculates summary statistics evaluating:
	- the model fits
	- the results of hierarchical partitioning
- *confidenceIntervals.R* --- calculates confidence intervals using the 
-   **plots** --- contains scripts for each of the plots used in the paper
    (except for figure 1, which was produced in Tikz)

    -   *figure2.R* --- Model output compared to observed data for weekly RSV-confirmed hospitalizations in children aged 0–2 in metropolitan Perth over the period 2000–2013
    -   *figure3.R* --- Effective (Reff) and basic (R0) reproduction numbers for the sinusoidal and weather models over time
    -   *figure4,5,sup.R* --- Model output of:
		- weekly RSV-confirmed hospitalizations in children aged 0–2 in metropolitan Perth over the period 2000–2020
		- weekly RSV-confirmed hospitalizations in children aged 0–2 in metropolitan Perth over the period 2020–mid 2021 under different lockdown scenarios
		- effective reproduction numbers over the period 2020–mid 2021 under different lockdown scenarios

[<h2 id="section.models"> models </h2>](#top)

- *fits.backwards.selection.Rdata* --- contains the fitted nested weather models and their AIC values (from **scripts**/*fit.BackwardsSelection.R*).
-   *fits.base.CV.Rdata* --- contains the 5 fitted base models corresponding to the five partitions used during cross-validation (from **scripts**/*fit.R)*
- *fit.base.Rdata* --- contains the base model fitted to the entire time series; used to produce Figures 3-5 (from **scripts**/*fit.R)*
-   *fits.weather.CV.Rdata* --- contains the 5 fitted weather models corresponding to the five partitions used during cross-validation (from
    **scripts**/*fit.R*)
- *fit.weather.Rdata* --- contains the weather model fitted to the entire time series; used to produce Figures 3-5 (from **scripts**/*fit.R)*
-   *diagnostics.fits.base.Rdata* --- contains different fits of the base model given different starting parameters (from *scripts**/*diagnostics.R*)
-   *diagnostics.fits.weather.Rdata* --- contains different fits of the weather model given different starting parameters (from *scripts**/*diagnostics.R*)

[<h2 id="section.plots"> plots </h2>](#top)

-   *figure2.tif* --- Model output compared to observed data for weekly RSV-confirmed hospitalizations in children aged 0–2 in metropolitan Perth over the period 2000–2013
-   *figure3.tif* --- Effective and basic reproduction numbers for the sinusoidal and weather models over time
-   *figure4.tif* --- Model output of weekly RSV-confirmed hospitalizations in children aged 0–2 in metropolitan Perth over the period 2000–2020
-   *figure5.tif* --- Model output of weekly RSV-confirmed hospitalizations in children aged 0–2 in metropolitan Perth over the period 2020–mid 2021
- *figure6.tif* --- S1 Figure. Effective reproduction numbers for the sinusoidal and weather models following COVID-19-associated NPI