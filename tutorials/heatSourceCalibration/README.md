# Heat-source calibration tutorial

This tutorial calibrates the projected depth-distribution closure used by
AdditiveFOAM heat sources. The supplied worked example uses a
`projectedGaussian` source with SS316L. The laser spot size (D4sigma) was
measured to be 109.69 microns.

## Installation and setup

Build AdditiveFOAM against OpenFOAM 14, then source both environments:

```bash
source /path/to/OpenFOAM-14/etc/bashrc
source /path/to/AdditiveFOAM/etc/bashrc
```

Copy the tutorial before running it:

```bash
mkdir -p "$FOAM_RUN/AdditiveFOAM"
cp -r "$ADDITIVEFOAM_TUTORIALS/heatSourceCalibration" \
    "$FOAM_RUN/AdditiveFOAM/heatSourceCalibration"
cd "$FOAM_RUN/AdditiveFOAM/heatSourceCalibration"
```

The repository installation creates the Python environment at
`$ADDITIVEFOAM_PROJECT_DIR/.venv`. It is activated by the AdditiveFOAM
environment:

```bash
calibrateHeatSource --help
```

Do not run a calibration inside `$ADDITIVEFOAM_TUTORIALS`; the campaign writes
cases and reports beneath the tutorial directory.

## Run the calibration

Review `system/decomposeParDict` in the template before starting. The supplied
research configuration requests 6 MPI ranks for each CFD case, evaluates ten
trial values for each of five experiments, and uses 2,000 posterior draws.

```bash
calibrateHeatSource --config config.yml
```

The configuration resolves relative paths from the location of `config.yml`.
Environment variables and `~` are supported in the three `paths` entries.

Generated output is contained in `campaign/`:

```text
campaign/
├── cases/
│   └── P187p5_V500_D109p69/
│       ├── B0/
│       ├── B4p5/
│       └── B9/
├── simulations.yml
├── calibration_state.yml
├── calibration_fit.yml
└── reports/
    ├── calibration_report.pdf
    └── calibration_summary.csv
```

Successful cases retain their rendered inputs, solver log, post-processing
output, and `metrics.yml`. With `keep_successful: false`, processor and numeric
time directories are removed after their melt-pool dimensions are recorded.

## Material-derived liquidus

The case includes:

```foam
#include "$ADDITIVEFOAM_ETC/materials/SS316L.cfg"
```

The configuration selects `melt_pool_isovalue: liquidus`. The calibration
command expands `constant/transportProperties` with `foamDictionary`, obtains
the temperature paired with `alpha.solid = 0`, and reads the matching file from
`postProcessing/meltPoolDimensions`. The SS316L liquidus is therefore not
duplicated in `config.yml`, `heatSourceDict`, or `controlDict`.

## Shared projected-source parameter

The template exposes the trial value once in `constant/heatSourceDict`:

```foam
calibrationA 0.0;
calibrationB <<B>>;
```

The template references these aliases directly from the projected Gaussian
coefficients:

```foam
projectedGaussianCoeffs
{
    dimensions (54.845e-6 54.845e-6 5.0e-5);
    A $calibrationA;
    B $calibrationB;
}
```

Power, speed, end time, and write interval are rendered independently.
