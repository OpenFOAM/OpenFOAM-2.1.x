/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "CarnahanStarlingRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(CarnahanStarling, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        CarnahanStarling,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::CarnahanStarling::CarnahanStarling
(
    const dictionary& dict
)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::CarnahanStarling::~CarnahanStarling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::CarnahanStarling::g0
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMax
) const
{

    return
        1.0/(1.0 - alpha1)
      + 3.0*alpha1/(2.0*sqr(1.0 - alpha1))
      + sqr(alpha1)/(2.0*pow(1.0 - alpha1, 3));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::CarnahanStarling::g0prime
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMax
) const
{
    return
        - alpha1/sqr(1.0 - alpha1)
        + (3.0*(1.0 - alpha1) + 6.0*sqr(alpha1))/(2.0*(1.0 - alpha1))
        + (2.0*alpha1*(1.0 - alpha1) + 3.0*pow(alpha1, 3))
         /(2.0*pow(1.0 - alpha1, 4));
}


// ************************************************************************* //
