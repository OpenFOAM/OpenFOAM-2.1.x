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

#include "mathematicalConstants.H"
#include "universalConstants.H"
#include "electromagneticConstants.H"
#include "atomicConstants.H"

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const char* const Foam::constant::atomic::group = "atomic";


const Foam::dimensionedScalar Foam::constant::atomic::alpha
(
    dimensionedConstant
    (
        group,
        "alpha",
        dimensionedScalar
        (
            "alpha",
            sqr(constant::electromagnetic::e)
           /(
                dimensionedScalar("C", dimless, 2.0)
               *constant::electromagnetic::epsilon0
               *constant::universal::h
               *constant::universal::c
            )
        )
    )
);


const Foam::dimensionedScalar Foam::constant::atomic::Rinf
(
    dimensionedConstant
    (
        group,
        "Rinf",
        dimensionedScalar
        (
            "Rinf",
            sqr(alpha)*me*constant::universal::c
           /(dimensionedScalar("C", dimless, 2.0)*constant::universal::h)
        )
    )
);


const Foam::dimensionedScalar Foam::constant::atomic::a0
(
    dimensionedConstant
    (
        group,
        "a0",
        dimensionedScalar
        (
            "a0",
            alpha
           /(
               dimensionedScalar("C", dimless, 4.0*constant::mathematical::pi)
              *Rinf
           )
        )
    )
);


const Foam::dimensionedScalar Foam::constant::atomic::re
(
    dimensionedConstant
    (
        group,
        "re",
        dimensionedScalar
        (
            "re",
            sqr(constant::electromagnetic::e)
           /(
                dimensionedScalar("C", dimless, 4.0*constant::mathematical::pi)
               *constant::electromagnetic::epsilon0
               *me
               *sqr(constant::universal::c)
            )
        )
    )
);


const Foam::dimensionedScalar Foam::constant::atomic::Eh
(
    dimensionedConstant
    (
        group,
        "Eh",
        dimensionedScalar
        (
            "Eh",
            dimensionedScalar("C", dimless, 2.0)
           *Rinf*constant::universal::h*constant::universal::c
        )
    )
);


// ************************************************************************* //
