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
    Foam::exponentialSolidThermo

Description
    Exponential thermodynamics package for solids templated into the rhoType.

SourceFiles
    exponentialSolidThermoI.H
    exponentialSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialSolidThermo_H
#define exponentialSolidThermo_H

#include "scalar.H"
#include "dictionary.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class rhoType> class exponentialSolidThermo;

template<class rhoType>
Ostream& operator<<
(
    Ostream&,
    const exponentialSolidThermo<rhoType>&
);


/*---------------------------------------------------------------------------*\
                         Class exponentialSolidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class rhoType>
class exponentialSolidThermo
:
    public rhoType
{

    // Private data

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar Hf_;


    // Private Member Functions

        //- Check given temperature is within the range of the fitted coeffs
        inline void checkT(const scalar T) const;

        //- Integrate Cp expression
        inline scalar integrateCp(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline exponentialSolidThermo
        (
            const rhoType& st,
            const scalar c0,
            const scalar n0,
            const scalar Tref,
            const scalar Hf
        );

        //- Construct from dictionary
        exponentialSolidThermo(const dictionary&);

        //- Construct as a named copy
        inline exponentialSolidThermo
        (
            const word&,
            const exponentialSolidThermo&
        );


    // Member Functions

        //- Heat capacity at constant pressure [J/(kg K)]
        inline scalar Cp(const scalar T) const;

        //- Enthalpy [J/kmol]
        inline scalar h(const scalar T) const;

        //- Sensible enthalpy [J/kg]
        inline scalar hs(const scalar T) const;

        //- Formation enthalpy [J/kg]
        inline scalar hf() const;


    // Ostream Operator

        friend Ostream& operator<< <rhoType>
        (
            Ostream&,
            const exponentialSolidThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "exponentialSolidThermoI.H"
#   include "exponentialSolidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
