if ((runTime.value() - endTime) > small)
{
    scalar maxDi =
        runTime.controlDict().lookupOrDefault<scalar>("maxDi", 10.0);

    const volScalarField::Internal DiNumvf
    (
        fvc::surfaceSum
        (
            mesh.magSf()
           *fvc::interpolate(kappa)
           *mesh.surfaceInterpolation::deltaCoeffs()
        )()()
       /(mesh.V()*rho.value()*Cp())
       *runTime.deltaT()
    );

    const scalar DiNum = gMax(DiNumvf);

    const scalar maxDeltaTFact = maxDi/(DiNum + small);
    
    const scalar deltaTFact =
        min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    const scalar deltaTDi =
        min(deltaTFact*runTime.deltaTValue(), maxDeltaT);

    runTime.setDeltaT(deltaTDi);
}
else
{
    runTime.setDeltaT(deltaT0);
}

Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
