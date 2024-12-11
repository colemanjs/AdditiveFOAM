forAll(mesh.cells(), celli)
{
    if (pos(sources.qDot()[celli] - small))
    {
        alpha1[celli] = 1.0;
    }

    const scalar& x = T[celli];
    const scalar& f = alpha1[celli];
          
    kappa[celli] = f*kappa1.value(x) + (1.0 - f)*kappa2.value(x);
    Cp[celli] = f*Cp1.value(x) + (1.0 - f)*Cp2.value(x);
}

alpha1.correctBoundaryConditions();
Cp.correctBoundaryConditions();
kappa.correctBoundaryConditions();
