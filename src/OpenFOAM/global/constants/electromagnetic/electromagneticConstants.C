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

const char* const Foam::constant::electromagnetic::group = "electromagnetic";


const Foam::dimensionedScalar Foam::constant::electromagnetic::mu0
(
    dimensionedConstant
    (
        group,
        "mu0",
        dimensionedScalar
        (
            "mu0",
            dimensionSet(1, 1, -2, 0, 0, -2, 0),
            4.0*constant::mathematical::pi*1e-07
        )
    )
);


const Foam::dimensionedScalar Foam::constant::electromagnetic::epsilon0
(
    dimensionedConstant
    (
        group,
        "epsilon0",
        dimensionedScalar
        (
            "epsilon0",
            dimensionedScalar("C", dimless, 1.0)
           /(mu0*sqr(constant::universal::c))
        )
    )
);


const Foam::dimensionedScalar Foam::constant::electromagnetic::Z0
(
    dimensionedConstant
    (
        group,
        "Z0",
        dimensionedScalar
        (
            "Z0",
            mu0*constant::universal::c
        )
    )
);


const Foam::dimensionedScalar Foam::constant::electromagnetic::kappa
(
    dimensionedConstant
    (
        group,
        "kappa",
        dimensionedScalar
        (
            "kappa",
            dimensionedScalar
            (
                "C",
                dimless,
                1.0/(4.0*constant::mathematical::pi)
            )
           /epsilon0
        )
    )
);


const Foam::dimensionedScalar Foam::constant::electromagnetic::G0
(
    dimensionedConstant
    (
        group,
        "G0",
        dimensionedScalar
        (
            "G0",
            dimensionedScalar("C", dimless, 2)*sqr(e)/constant::universal::h
        )
    )
);


const Foam::dimensionedScalar Foam::constant::electromagnetic::KJ
(
    dimensionedConstant
    (
        group,
        "KJ",
        dimensionedScalar
        (
            "KJ",
            dimensionedScalar("C", dimless, 2)*e/constant::universal::h
        )
    )
);


const Foam::dimensionedScalar Foam::constant::electromagnetic::phi0
(
    dimensionedConstant
    (
        group,
        "phi0",
        dimensionedScalar
        (
            "phi0",
            constant::universal::h/(dimensionedScalar("C", dimless, 2)*e)
        )
    )
);


const Foam::dimensionedScalar Foam::constant::electromagnetic::RK
(
    dimensionedConstant
    (
        group,
        "RK",
        dimensionedScalar
        (
            "RK",
            constant::universal::h/sqr(e)
        )
    )
);


// ************************************************************************* //
