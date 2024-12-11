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
    Foam::heatSourceModels::Gaussian

Description
    Gaussian heat source model for additive manufacting

SourceFiles
    Gaussian.C

\*---------------------------------------------------------------------------*/

#ifndef Gaussian_H
#define Gaussian_H

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{

/*---------------------------------------------------------------------------*\
                       Class Gaussian
\*---------------------------------------------------------------------------*/

class Gaussian
:
    public heatSourceModel
{
    // Private Data
    
        //- Pointer to mesh information
        const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("Gaussian");


    // Constructors

        //- Construct from components
        Gaussian
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Gaussian()
    {}


    // Member Functions

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatSourceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
