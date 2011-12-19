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

#include "thermalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMedia
{
   defineTypeNameAndDebug(thermalModel, 0);
   defineRunTimeSelectionTable(thermalModel, pZone);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMedia::thermalModel::thermalModel(const porousZone& pZone)
:
    pZone_(pZone),
    thermalCoeffs_(pZone.dict().subDictPtr("thermalModel"))
{}


Foam::porousMedia::thermalModel::thermalModel
(
    const porousZone& pZone,
    const dictionary& thermalCoeffs
)
:
    pZone_(pZone),
    thermalCoeffs_(thermalCoeffs)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::porousMedia::thermalModel::~thermalModel()
{}


// ************************************************************************* //
