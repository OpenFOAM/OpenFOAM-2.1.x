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
    Foam::constRho

Description
    rho constant

SourceFiles
    constRhoI.H
    constRho.C

\*---------------------------------------------------------------------------*/

#ifndef constRho_H
#define constRho_H

#include "dictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constRho Declaration
\*---------------------------------------------------------------------------*/

class constRho
{
    // Private data

        //- Density [kg/m3]
        scalar rho_;


public:

    // Constructors

        //- Construct from components
        inline constRho(scalar rho);

        //- Construct from dictionary
        //constRho(Istream&);
        constRho(const dictionary& dict);

        //- Construct as named copy
        inline constRho(const constRho&);

        //- Construct and return a clone
        inline autoPtr<constRho> clone() const;

        // Selector from dictionary
        inline static autoPtr<constRho> New(const dictionary& is);


    // Member functions

        //- Return density [kg/m3]
        inline scalar rho(scalar T) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const constRho&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constRhoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
