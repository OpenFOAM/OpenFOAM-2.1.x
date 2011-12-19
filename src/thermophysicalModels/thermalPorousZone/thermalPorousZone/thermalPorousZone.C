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

\*----------------------------------------------------------------------------*/

#include "thermalPorousZone.H"
#include "basicThermo.H"
#include "volFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPorousZone::thermalPorousZone
(
    const keyType& key,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    porousZone(key, mesh, dict),
    model_(porousMedia::thermalModel::New(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermalPorousZone::addEnthalpySource
(
    const basicThermo& thermo,
    const volScalarField& rho,
    fvScalarMatrix& hEqn
) const
{
    if (model_.valid())
    {
        model_->addEnthalpySource(thermo, rho, hEqn);
    }
}


// ************************************************************************* //
