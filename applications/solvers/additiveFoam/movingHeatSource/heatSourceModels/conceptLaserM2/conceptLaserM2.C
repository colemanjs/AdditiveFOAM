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

#include "conceptLaserM2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(conceptLaserM2, 0);
    addToRunTimeSelectionTable(heatSourceModel, conceptLaserM2, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::conceptLaserM2::conceptLaserM2
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh)
{
    // NOTE: Hardcoded adjustment for ellipsticity (about 10%)
    dimensions_ =
        vector(0.95 * dimensions_.x(), 0.95 * dimensions_.y(), dimensions_.z());
    staticDimensions_ = dimensions_;

    sigma0 = 0.5 * dimensions_.x();    
    d0 = dimensions_.z();
    
    Info << sigma0 << tab << d0 << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::conceptLaserM2::weight(const vector& x)
{
    scalar z = Foam::mag(x.z());

    scalar r = std::sqrt(x.x() * x.x() + x.y() * x.y());
        
    scalar x0 = 2.0 * std::exp(-0.5 * Foam::sqr(r / sigma0));
    
    scalar s0 = std::exp(-3.0 * std::pow(mag(z/d0), n0));
           
    return x0 * s0;
}

inline Foam::scalar
Foam::heatSourceModels::conceptLaserM2::ai(scalar s, scalar d, scalar n)
{
    const scalar t1 =
        4.0 * pi * s * s * d * Foam::tgamma(1.0 / n)
      / ( n * std::pow(3.0, 1.0 / n) );
      
    return t1;
}


inline Foam::dimensionedScalar
Foam::heatSourceModels::conceptLaserM2::V0()
{
    // NOTE: THE CALIBRATED SHAPE FUNCTION IS HARDCODED FOR FOR NOW
    d0 = dimensions_.z();

    const scalar x = max(d0 / staticDimensions_.x(), 1.0);
        
    n0 = min(max(4.2*std::log2(x), 1.0), 9.0);

    Info << "x: " << x << tab << "n: " << n0 << endl;
                
    n0 = std::pow(2.0, n0);
    
    a0 = ai(sigma0, d0, n0);
    
    return dimensionedScalar("V0", dimVolume, a0);
}



bool Foam::heatSourceModels::conceptLaserM2::read()
{
    if (heatSourceModel::read())
    {
        //- Mandatory entries        
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
