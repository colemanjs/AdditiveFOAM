Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField alpha1
(
    IOobject
    (
        IOobject::groupName("alpha", "solid"),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

movingHeatSourceModel sources(mesh);

#include "readTransportProperties.H"

volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cp", dimensionSet(0, 2, -2, -1, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("kappa", dimensionSet(1, 1, -3, -1, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

#include "updateProperties.H"
