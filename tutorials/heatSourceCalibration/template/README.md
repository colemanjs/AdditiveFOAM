# Calibration case template

This is the OpenFOAM-14 AdditiveFOAM case copied for every trial value in the
heat-source calibration campaign. It uses a `projectedGaussian` heat source
with SS316L. The laser spot size (D4sigma) was measured to be 109.69 microns.

The calibration command expects this full case structure:

```text
0/
constant/
system/
Allrun
Allclean
```

The required renderer placeholders are:

- `constant/heatSourceDict`: `<<B>>`, exposed once as `calibrationB`
- `constant/scanPath`: `<<power>>`, `<<velocity>>`
- `system/controlDict`: `<<endTime>>`, `<<writeInterval>>`

OpenFOAM dictionary aliases propagate `$calibrationB` into the projected
Gaussian closure. A template for another projected source can use the same
alias for every applicable `B` coefficient.

The SS316L material configuration supplies `thermoPath`, emissivity, and the
Marangoni coefficient. Transient source depth and `meltPoolDimensions` obtain
their default isovalues from that material configuration.
