/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::pressureGradientExplicitSource

Description
    Creates a pressure gradient source

    Note: Currently only handles kinematic pressure

    Sources described by:

        pressureGradientExplicitSourceCoeffs
        {
            UName       U;          // name of velocity field
            Ubar        (10.0 0 0); // desired average velocity
            gradPini    gradPini [0 2 -2 0 0] 0; // initial pressure gradient
            flowDir     (1 0 0);    // flow direction
        }
        

SourceFiles
    pressureGradientExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef pressureGradientExplicitSource_H
#define pressureGradientExplicitSource_H

#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "basicSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pressureGradientExplicitSource Declaration
\*---------------------------------------------------------------------------*/

class pressureGradientExplicitSource
:
    public basicSource
{
    // Private data

        //- Average velocity
        vector Ubar_;

        //- Initial pressure gradient
        dimensionedScalar gradPini_;

        //- Pressure gradient
        dimensionedScalar gradP_;

        //- Flow direction
        vector flowDir_;


    // Private Member Functions

        //- Write the pressure gradient to file (for restarts etc)
        void writeGradP() const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);

        //- Disallow default bitwise copy construct
        pressureGradientExplicitSource(const pressureGradientExplicitSource&);

        //- Disallow default bitwise assignment
        void operator=(const pressureGradientExplicitSource&);


public:

    //- Runtime type information
    TypeName("pressureGradientExplicitSource");


    // Constructors

        //- Construct from explicit source name and mesh
        pressureGradientExplicitSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to equation
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
