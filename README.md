# Summary
This package implements several methods for epidemiological forecasting.
It also implements recalibration methods for epidemiological forecasting, described in this [paper](https://arxiv.org/abs/2112.06305).

# Recalibration
The recalibration methods are implemented in `epiforecast/R/calibrate.R` and scripts demonstrating their use on the FluSight Ensemble forecasters are located in `calibration-experiments`.

To run the recalibration experiments, perform the following commands in the `calibration-experiments` directory:
1. Install the FluSight repo:
`git clone https://github.com/FluSightNetwork/cdc-flusight-ensemble.git`
1. Create the component forecast file: `Rscript save_historical_forecasts.R`
1. Run desired `calibrate_*.R` scripts:
    1. `calibrate_cv.R` recalibrates each of the FluSight Network component forecasters in an out-of-sample fashion.
    1. `calibrate_ensemble.R` recalibrates the FluSight Network ensemble itself out-of-sample (see section 3.3).
    1. `calibrate_sensitivity.R` recalibrates the FluSight Network component forecasters and experiments with different training window sizes (see sections 2.6, 3.1).
