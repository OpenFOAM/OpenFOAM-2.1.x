/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::pressurePerfectGas

Description
    Perfect gas equation of state usign a reference pressure

SourceFiles
    pressurePerfectGasI.H
    pressurePerfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef pressurePerfectGas_H
#define pressurePerfectGas_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pressurePerfectGas Declaration
\*---------------------------------------------------------------------------*/

class pressurePerfectGas
:
    public specie
{
    // Private data

        //- Reference pressure
        scalar pRef_;


public:

    // Constructors

        //- Construct from components
        inline pressurePerfectGas(const specie& sp);

        //- Construct from Istream
        pressurePerfectGas(Istream&);

        //- Construct from dictionary
        pressurePerfectGas(const dictionary& dict);

        //- Construct as named copy
        inline pressurePerfectGas(const word& name, const pressurePerfectGas&);

        //- Construct and return a clone
        inline autoPtr<pressurePerfectGas> clone() const;

        // Selector from Istream
        inline static autoPtr<pressurePerfectGas> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<pressurePerfectGas> New(const dictionary& dict);


    // Member functions

        // Fundamental properties

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const pressurePerfectGas&);
        inline void operator-=(const pressurePerfectGas&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend pressurePerfectGas operator+
        (
            const pressurePerfectGas&,
            const pressurePerfectGas&
        );

        inline friend pressurePerfectGas operator-
        (
            const pressurePerfectGas&,
            const pressurePerfectGas&
        );

        inline friend pressurePerfectGas operator*
        (
            const scalar s,
            const pressurePerfectGas&
        );

        inline friend pressurePerfectGas operator==
        (
            const pressurePerfectGas&,
            const pressurePerfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const pressurePerfectGas&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pressurePerfectGasI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
