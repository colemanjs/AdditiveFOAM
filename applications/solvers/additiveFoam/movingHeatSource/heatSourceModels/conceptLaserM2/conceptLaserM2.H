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

Class
    Foam::heatSourceModels::conceptLaserM2

Description
    conceptLaserM2 heat source model for additive manufacting

SourceFiles
    conceptLaserM2.C

\*---------------------------------------------------------------------------*/

#ifndef conceptLaserM2_H
#define conceptLaserM2_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                       Class conceptLaserM2
\*---------------------------------------------------------------------------*/

class conceptLaserM2
:
    public heatSourceModel
{
    // Private Data
    
        //- Pointer to mesh information
        const fvMesh& mesh_;
        
        //- Shape coefficient of conceptLaserM2 distribution        
        scalar d0;
        
        scalar sigma0;
                
        scalar a0;
        
        scalar n0;

public:

    //- Runtime type information
    TypeName("conceptLaserM2");


    // Constructors

        //- Construct from components
        conceptLaserM2
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~conceptLaserM2()
    {}


    // Member Functions

        inline virtual scalar weight(const vector& d);


        inline scalar ai(scalar s, scalar d, scalar m);
        
        inline virtual dimensionedScalar V0();

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
