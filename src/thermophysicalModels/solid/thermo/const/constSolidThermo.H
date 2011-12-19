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
    Foam::constSolidThermo

Description
    Constant properties thermodynamics package for solids


SourceFiles
    constSolidThermoI.H
    constSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constSolidThermo_H
#define constSolidThermo_H

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class rhoType> class constSolidThermo;

template<class rhoType>
Ostream& operator<<
(
    Ostream&,
    const constSolidThermo<rhoType>&
);


/*---------------------------------------------------------------------------*\
                      Class constSolidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class rhoType>
class constSolidThermo
:
    public rhoType
{
    // Private data

        scalar Cp_;
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline constSolidThermo
        (
            const rhoType& st,
            const scalar Cp,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        constSolidThermo(const dictionary&);

        //- Construct as named copy
        inline constSolidThermo(const word&, const constSolidThermo&);

        //- Construct and return a clone
        inline autoPtr<constSolidThermo> clone() const;

        //- Selector from Istream
        inline static autoPtr<constSolidThermo> New(dictionary& is);


    // Member Functions

        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar T) const;

            //- Enthalpy [J/kg]
            inline scalar h(const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar T) const;

            //- Enthalpy of formation [J/kg]
            inline scalar hf() const;


    // IOstream Operators

        friend Ostream& operator<< <rhoType>
        (
            Ostream&,
            const constSolidThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "constSolidThermoI.H"
#   include "constSolidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
