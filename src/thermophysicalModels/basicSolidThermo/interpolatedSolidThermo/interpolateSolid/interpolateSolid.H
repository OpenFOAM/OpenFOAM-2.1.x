/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::interpolateSolid

Description
    Helping class for T-interpolated solid thermo.

SourceFiles
    interpolateSolid.C

\*---------------------------------------------------------------------------*/

#ifndef interpolateSolid_H
#define interpolateSolid_H
#include "volFields.H"
#include "dictionary.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class interpolateSolid Declaration
\*---------------------------------------------------------------------------*/

class interpolateSolid
{
protected:

    // Protected data

        Field<scalar> TValues_;

        Field<scalar> rhoValues_;

        Field<scalar> cpValues_;

        Field<scalar> HfValues_;

        Field<scalar> emissivityValues_;

        Field<scalar> kappaValues_;

        Field<scalar> sigmaSValues_;


public:

    // Constructors

        //- Construct from mesh
        interpolateSolid(const dictionary&);


    //- Destructor
    virtual ~interpolateSolid();


    // Member Functions

        // I-O

            //- Read the interpolateSolid properties
            bool read(const dictionary& dict);

            //- Write the interpolateSolid properties
            bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
