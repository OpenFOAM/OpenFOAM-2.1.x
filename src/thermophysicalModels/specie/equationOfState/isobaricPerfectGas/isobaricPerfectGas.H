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
    Foam::isobaricPerfectGas

Description
    Perfect gas equation of state using a reference pressure
    rather than the local pressure.

SourceFiles
    isobaricPerfectGasI.H
    isobaricPerfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef isobaricPerfectGas_H
#define isobaricPerfectGas_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class isobaricPerfectGas Declaration
\*---------------------------------------------------------------------------*/

class isobaricPerfectGas
:
    public specie
{
    // Private data

        //- Reference pressure
        scalar pRef_;


public:

    // Constructors

        //- Construct from components
        inline isobaricPerfectGas(const specie& sp);

        //- Construct from Istream
        isobaricPerfectGas(Istream&);

        //- Construct from dictionary
        isobaricPerfectGas(const dictionary& dict);

        //- Construct as named copy
        inline isobaricPerfectGas(const word& name, const isobaricPerfectGas&);

        //- Construct and return a clone
        inline autoPtr<isobaricPerfectGas> clone() const;

        // Selector from Istream
        inline static autoPtr<isobaricPerfectGas> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<isobaricPerfectGas> New(const dictionary& dict);


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

        inline void operator+=(const isobaricPerfectGas&);
        inline void operator-=(const isobaricPerfectGas&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend isobaricPerfectGas operator+
        (
            const isobaricPerfectGas&,
            const isobaricPerfectGas&
        );

        inline friend isobaricPerfectGas operator-
        (
            const isobaricPerfectGas&,
            const isobaricPerfectGas&
        );

        inline friend isobaricPerfectGas operator*
        (
            const scalar s,
            const isobaricPerfectGas&
        );

        inline friend isobaricPerfectGas operator==
        (
            const isobaricPerfectGas&,
            const isobaricPerfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const isobaricPerfectGas&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "isobaricPerfectGasI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
