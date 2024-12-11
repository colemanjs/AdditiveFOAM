IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const int PolySize(3);
const dictionary& solidDict(transportProperties.subDict("solid"));
Polynomial<PolySize> kappa1(solidDict.lookup("kappa"));
Polynomial<PolySize> Cp1(solidDict.lookup("Cp"));

const dictionary& powderDict(transportProperties.subDict("powder"));
Polynomial<PolySize> kappa2(powderDict.lookup("kappa"));
Polynomial<PolySize> Cp2(powderDict.lookup("Cp"));

const dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);
