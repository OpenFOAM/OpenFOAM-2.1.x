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

Description
    Fundamental dimensioned constants

\*---------------------------------------------------------------------------*/

#include "fundamentalConstants.H"

#include "universalConstants.H"
#include "electromagneticConstants.H"
#include "atomicConstants.H"
#include "physicoChemicalConstants.H"

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Universal constants

const Foam::dimensionedScalar Foam::constant::universal::c
(
    dimensionedConstant(universal::group, "c")
);


const Foam::dimensionedScalar Foam::constant::universal::G
(
    dimensionedConstant(universal::group, "G")
);


const Foam::dimensionedScalar Foam::constant::universal::h
(
    dimensionedConstant(universal::group, "h")
);


// Electromagnetic

const Foam::dimensionedScalar Foam::constant::electromagnetic::e
(
    dimensionedConstant(electromagnetic::group, "e")
);


// Atomic

const Foam::dimensionedScalar Foam::constant::atomic::me
(
    dimensionedConstant(atomic::group, "me")
);


const Foam::dimensionedScalar Foam::constant::atomic::mp
(
    dimensionedConstant(atomic::group, "mp")
);


// Physico-chemical

const Foam::dimensionedScalar Foam::constant::physicoChemical::mu
(
    dimensionedConstant(physicoChemical::group, "mu")
);


const Foam::dimensionedScalar Foam::constant::physicoChemical::NA
(
//    dimensionedConstant(physicoChemical::group, "NA")
    dimensionedConstant
    (
        physicoChemical::group,
        "NA",
        dimensionedScalar
        (
            "NA",
            dimless/dimMoles,
            6.0221417930e+23
        )
    )
);


const Foam::dimensionedScalar Foam::constant::physicoChemical::k
(
    dimensionedConstant(physicoChemical::group, "k")
);


// Standard

const Foam::dimensionedScalar Foam::constant::standard::Pstd
(
    dimensionedConstant("standard", "Pstd")
);


const Foam::dimensionedScalar Foam::constant::standard::Tstd
(
    dimensionedConstant("standard", "Tstd")
);


// ************************************************************************* //
