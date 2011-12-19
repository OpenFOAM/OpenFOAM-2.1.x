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

\*---------------------------------------------------------------------------*/

#include "perfectGas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::perfectGas::perfectGas(const specie& sp)
:
    specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::perfectGas::perfectGas(const word& name, const perfectGas& pg)
:
    specie(name, pg)
{}


inline Foam::autoPtr<Foam::perfectGas> Foam::perfectGas::clone() const
{
    return autoPtr<perfectGas>(new perfectGas(*this));
}


inline Foam::autoPtr<Foam::perfectGas> Foam::perfectGas::New(Istream& is)
{
    return autoPtr<perfectGas>(new perfectGas(is));
}


inline Foam::autoPtr<Foam::perfectGas> Foam::perfectGas::New
(
    const dictionary& dict
)
{
    return autoPtr<perfectGas>(new perfectGas(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::perfectGas::rho(scalar p, scalar T) const
{
    return p/(R()*T);
}


inline Foam::scalar Foam::perfectGas::psi(scalar, scalar T) const
{
    return 1.0/(R()*T);
}


inline Foam::scalar Foam::perfectGas::Z(scalar, scalar) const
{
    return 1.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::perfectGas::operator+=(const perfectGas& pg)
{
    specie::operator+=(pg);
}


inline void Foam::perfectGas::operator-=(const perfectGas& pg)
{
    specie::operator-=(pg);
}


inline void Foam::perfectGas::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::perfectGas Foam::operator+
(
    const perfectGas& pg1,
    const perfectGas& pg2
)
{
    return perfectGas
    (
        static_cast<const specie&>(pg1)
      + static_cast<const specie&>(pg2)
    );
}


inline Foam::perfectGas Foam::operator-
(
    const perfectGas& pg1,
    const perfectGas& pg2
)
{
    return perfectGas
    (
        static_cast<const specie&>(pg1)
      - static_cast<const specie&>(pg2)
    );
}


inline Foam::perfectGas Foam::operator*
(
    const scalar s,
    const perfectGas& pg
)
{
    return perfectGas(s*static_cast<const specie&>(pg));
}


inline Foam::perfectGas Foam::operator==
(
    const perfectGas& pg1,
    const perfectGas& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
