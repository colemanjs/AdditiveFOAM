/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "heatSourceModel.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatSourceModel, 0);
    defineRunTimeSelectionTable(heatSourceModel, dictionary);
}

const Foam::word Foam::heatSourceModel::heatSourceDictName
(
    "heatSourceDict"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::heatSourceModel::createIOobject
(
    const dictionary& dict,
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        dict.name(),
        mesh.time().constant(),
        mesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModel::heatSourceModel
(
    const word& type,
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(dict, mesh)),

    sourceName_(sourceName),
    heatSourceDict_(dict),
    sourceDict_(heatSourceDict_.optionalSubDict(sourceName_)),
    heatSourceModelCoeffs_(sourceDict_.optionalSubDict(type + "Coeffs")),
    
    mesh_(mesh),
    absorptionModel_(nullptr),
    movingBeam_(nullptr)
{
    absorptionModel_ =
        absorptionModel::New(sourceName_, heatSourceDict_, mesh_);

    movingBeam_ =
        movingBeam::New(sourceName_, heatSourceDict_, mesh_.time());

    dimensions_ =
        heatSourceModelCoeffs_.lookup<vector>("dimensions");

    staticDimensions_ = dimensions_;

    transient_ =
        heatSourceModelCoeffs_.lookupOrDefault<Switch>("transient", false);

    isoValue_ =
        heatSourceModelCoeffs_.lookupOrDefault<scalar>("isoValue", great);

    nPoints_ = 
        heatSourceModelCoeffs_.lookupOrDefault<labelVector>
        (
            "nPoints",
            vector::one
        );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatSourceModel::qDot()
{
    tmp<volScalarField> tqDot
    (
        new volScalarField
        (
            IOobject
            (
                "qDot_",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Zero", dimPower/dimVolume, 0.0)
        )
    );
    volScalarField& qDot_ = tqDot.ref();

    // sample heat source distribution at desired resolution
    const scalar power_ = movingBeam_->power();

    if (power_ > small)
    {    
        const vector position_ = movingBeam_->position();

        // udpate the absorbed power and heat source normalization term
        const scalar aspectRatio = 
            dimensions_.z() / min(dimensions_.x(), dimensions_.y());

        dimensionedScalar absorbedPower
        (
            "etaP",
            dimPower,
            absorptionModel_->eta(aspectRatio)*power_
        );

        dimensionedScalar volume = V0();

        // integrate the heat source in each overlapping cell
        volScalarField weights
        (
            IOobject
            (
                "weights",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Zero", dimless, 0.0)          
        );

        const pointField& points = mesh_.points();
        
        treeBoundBox beamBb
        (
            position_ - 2.0*dimensions_,
            position_ + 2.0*dimensions_
        );

        forAll(mesh_.cells(), celli)
        {
            treeBoundBox cellBb(point::max, point::min);

            const labelList& vertices = mesh_.cellPoints()[celli];

            forAll(vertices, i)
            {
                cellBb.min() = min(cellBb.min(), points[vertices[i]]);
                cellBb.max() = max(cellBb.max(), points[vertices[i]]);
            }

            if (cellBb.overlaps(beamBb))
            {
                vector dx_ = cmptDivide(dimensions_, vector(nPoints_));
              
                labelVector nCellPoints =
                    max
                    (
                        cmptDivide(cellBb.span() + small*vector::one, dx_),
                        vector::one
                    );

                dx_ = cmptDivide(cellBb.span(), vector(nCellPoints));

                label pI = 0;
                
                List<point> d(cmptProduct(nCellPoints));

                for (label k=0; k < nCellPoints.z(); ++k)
                {
                    for (label j=0; j < nCellPoints.y(); ++j)
                    {
                        for (label i=0; i < nCellPoints.x(); ++i)
                        {
                            const point pt
                            (
                                cellBb.max()
                              - cmptMultiply
                                (
                                    vector(i + 0.5, j + 0.5, k + 0.5),
                                    dx_
                                )
                            );

                            treeBoundBox ptBb(pt - 0.5*dx_, pt + 0.5*dx_);

                            if (beamBb.overlaps(ptBb))
                            {
                                d[pI] = cmptMag(pt - position_);
                                pI++;
                            }                                
                        }
                    }
                }
                d.resize(pI);
                
                scalar dVi = cmptProduct(dx_);
                weights[celli] = weight(d) * dVi / mesh_.V()[celli];
            }
        }

        // stabilize numerical integration errors within 95% of applied power
        dimensionedScalar sumWeights = fvc::domainIntegrate(weights);
        scalar residual = (sumWeights / volume).value();

        if (mag(1 - residual) < 0.05)
        {
            volume = sumWeights;
        }

        qDot_ = absorbedPower * weights / volume;
    }

    return tqDot;
}


bool Foam::heatSourceModel::read()
{
    if (regIOobject::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
