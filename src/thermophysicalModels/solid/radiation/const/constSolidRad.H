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
    Foam::constSolidRad

Description
    Constant radiative properties

SourceFiles
    constSolidRadI.H
    constSolidRad.C

\*---------------------------------------------------------------------------*/

#ifndef constSolidRad_H
#define constSolidRad_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class thermo> class constSolidRad;

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const constSolidRad<thermo>&
);

/*---------------------------------------------------------------------------*\
                      Class constSolidRad Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class constSolidRad
:
    public  thermo
{
    // Private data

        //- Absorption coefficient
        scalar kappa_;

        //- Scattering coefficient
        scalar sigmaS_;

        //- Emissivity coefficient
        scalar emissivity_;


         //- Construct from components
        inline constSolidRad
        (
            const thermo& t,
            const scalar kappa,
            const scalar sigmaS,
            const scalar emissivity
        );


public:

    // Constructors

        //- Construct from Istream
        constSolidRad(const dictionary&);

        //- Construct as named copy
        inline constSolidRad(const word& name, const constSolidRad&);

        //- Construct and return a clone
        inline autoPtr<constSolidRad> clone() const;

        // Selector from Istream
        inline static autoPtr<constSolidRad> New(Istream& is);


    // Member functions

        //- Return absorption coefficient [1/m]
        inline scalar kappa(scalar T) const;

        //- Return scattering [1/m]
        inline scalar sigmaS(scalar T) const;

        //- Return  emissivity[]
        inline scalar emissivity(scalar T) const;


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const constSolidRad&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constSolidRadI.H"

#ifdef NoRepository
#   include "constSolidRad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
