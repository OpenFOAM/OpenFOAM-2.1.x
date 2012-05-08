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

\*---------------------------------------------------------------------------*/

#include "pressurePerfectGas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::pressurePerfectGas::pressurePerfectGas(const specie& sp)
:
    specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pressurePerfectGas::pressurePerfectGas
(
    const word& name,
    const pressurePerfectGas& pg
)
:
    specie(name, pg)
{}


inline Foam::autoPtr<Foam::pressurePerfectGas> Foam::pressurePerfectGas::
clone() const
{
    return autoPtr<pressurePerfectGas>(new pressurePerfectGas(*this));
}


inline Foam::autoPtr<Foam::pressurePerfectGas> Foam::pressurePerfectGas::New
(
    Istream& is
)
{
    return autoPtr<pressurePerfectGas>(new pressurePerfectGas(is));
}


inline Foam::autoPtr<Foam::pressurePerfectGas> Foam::pressurePerfectGas::New
(
    const dictionary& dict
)
{
    return autoPtr<pressurePerfectGas>(new pressurePerfectGas(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::pressurePerfectGas::rho(scalar p, scalar T) const
{
    return pRef_/(R()*T);
}


inline Foam::scalar Foam::pressurePerfectGas::psi(scalar, scalar T) const
{
    return 0.0;
}


inline Foam::scalar Foam::pressurePerfectGas::Z(scalar, scalar) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::pressurePerfectGas::operator+=(const pressurePerfectGas& pg)
{
    specie::operator+=(pg);
}


inline void Foam::pressurePerfectGas::operator-=(const pressurePerfectGas& pg)
{
    specie::operator-=(pg);
}


inline void Foam::pressurePerfectGas::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::pressurePerfectGas Foam::operator+
(
    const pressurePerfectGas& pg1,
    const pressurePerfectGas& pg2
)
{
    return pressurePerfectGas
    (
        static_cast<const specie&>(pg1)
      + static_cast<const specie&>(pg2)
    );
}


inline Foam::pressurePerfectGas Foam::operator-
(
    const pressurePerfectGas& pg1,
    const pressurePerfectGas& pg2
)
{
    return pressurePerfectGas
    (
        static_cast<const specie&>(pg1)
      - static_cast<const specie&>(pg2)
    );
}


inline Foam::pressurePerfectGas Foam::operator*
(
    const scalar s,
    const pressurePerfectGas& pg
)
{
    return pressurePerfectGas(s*static_cast<const specie&>(pg));
}


inline Foam::pressurePerfectGas Foam::operator==
(
    const pressurePerfectGas& pg1,
    const pressurePerfectGas& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
