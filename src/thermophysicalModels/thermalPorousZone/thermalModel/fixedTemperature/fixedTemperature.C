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

#include "fixedTemperature.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"
#include "volFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMedia
{
    defineTypeNameAndDebug(fixedTemperature, 0);

    addToRunTimeSelectionTable
    (
        thermalModel,
        fixedTemperature,
        pZone
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMedia::fixedTemperature::fixedTemperature(const porousZone& pZone)
:
    thermalModel(pZone),
    T_(readScalar(thermalCoeffs_.lookup("T")))
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::porousMedia::fixedTemperature::~fixedTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousMedia::fixedTemperature::addEnthalpySource
(
    const basicThermo& thermo,
    const volScalarField& rho,
    fvScalarMatrix& hEqn
) const
{
    const labelList& zones = pZone_.zoneIds();
    if (zones.empty() || T_ < 0.0)
    {
        return;
    }

    const fvMesh& mesh = pZone_.mesh();
    const scalarField& V = mesh.V();
    scalarField& hDiag = hEqn.diag();
    scalarField& hSource = hEqn.source();

    const scalarField T(hDiag.size(), T_);

    const scalar rate = 1e6;

    forAll(zones, zoneI)
    {
        const labelList& cells = mesh.cellZones()[zones[zoneI]];
        tmp<scalarField> h = thermo.h(T, cells);

        forAll(cells, i)
        {
            hDiag[cells[i]] += rate*V[cells[i]]*rho[cells[i]];
            hSource[cells[i]] +=
                rate*V[cells[i]]*rho[cells[i]]*h()[cells[i]];
        }
    }
}


// ************************************************************************* //
