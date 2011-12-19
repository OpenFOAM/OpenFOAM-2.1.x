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

#include "constRho.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::constRho::constRho
(
    scalar rho
)
:
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::constRho::constRho(const constRho& pg)
:
    rho_(pg.rho_)
{}


inline Foam::autoPtr<Foam::constRho> Foam::constRho::clone() const
{
    return autoPtr<constRho>(new constRho(*this));
}


inline Foam::autoPtr<Foam::constRho> Foam::constRho::New(const dictionary& is)
{
    return autoPtr<constRho>(new constRho(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::constRho::rho(scalar) const
{
    return rho_;
}


// ************************************************************************* //
