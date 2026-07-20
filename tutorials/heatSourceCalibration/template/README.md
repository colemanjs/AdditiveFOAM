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

- `constant/heatSourceDict`: `<<B>>` and `<<spotSize2Sigma>>`
- `constant/scanPath`: `<<power>>`, `<<velocity>>`
- `system/controlDict`: `<<endTime>>`, `<<writeInterval>>`

The projected Gaussian closure reads `<<B>>` directly. The
`$spotSize2Sigma` alias supplies both lateral dimensions from one rendered
value.

The calibration command converts `Spot_Size_microns` from the measured D4sigma
diameter to a 2sigma radius in metres and assigns it to both lateral source
dimensions.

The SS316L material configuration supplies `thermoPath`, emissivity, and the
Marangoni coefficient. Transient source depth and `meltPoolDimensions` obtain
their default isovalues from that material configuration.
