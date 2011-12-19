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
    Foam::exponentialSolidTransport

Description
    Exponential properties for solid heat transport
    Templated into a given thermodynamics package.

SourceFiles
    exponentialSolidTransportI.H
    exponentialSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialSolidTransport_H
#define exponentialSolidTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class thermo> class exponentialSolidTransport;

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const exponentialSolidTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                        Class exponentialSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class exponentialSolidTransport
:
    public thermo
{
    // Private data

        //- Constant thermal coefficient.
        scalar K0_;

        //- Exponent coefficient
        scalar n0_;

        //- Reference temperature
        scalar Tref_;


    // Private Member Functions

        //- Construct from components
        inline exponentialSolidTransport
        (
            const thermo& t,
            const scalar K0,
            const scalar n0,
            const scalar Tref
        );


public:

    // Constructors

        //- Construct as named copy
        inline exponentialSolidTransport
        (
            const word&,
            const exponentialSolidTransport&
        );


        //- Construct from dictionary
        exponentialSolidTransport(const dictionary&);


    // Member functions

        //- Thermal conductivity [W/mK]
        inline scalar K(const scalar T) const;

        //- Thermal diffusivity [kg/ms]
        inline scalar alpha(const scalar T) const;


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const exponentialSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "exponentialSolidTransportI.H"

#ifdef NoRepository
#   include "exponentialSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
