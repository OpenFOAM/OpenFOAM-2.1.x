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

#include "reactionTypes.H"
#include "makeReactionThermo.H"

#include "ArrheniusReactionRate.H"
#include "infiniteReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "FallOffReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(gasReaction, 0);
defineTemplateRunTimeSelectionTable(gasReaction, Istream);
defineTemplateRunTimeSelectionTable(gasReaction, dictionary);


// * * * * * * * * * * * * * Make CHEMKIN reactions  * * * * * * * * * * * * //

makeIRNReactions(gasThermoPhysics, ArrheniusReactionRate)
makeIRNReactions(gasThermoPhysics, infiniteReactionRate)
makeIRNReactions(gasThermoPhysics, LandauTellerReactionRate)
makeIRNReactions(gasThermoPhysics, thirdBodyArrheniusReactionRate)

makeIRReactions(gasThermoPhysics, JanevReactionRate)
makeIRReactions(gasThermoPhysics, powerSeriesReactionRate)

makePressureDependentReactions
(
    gasThermoPhysics,
    ArrheniusReactionRate,
    LindemannFallOffFunction
)

makePressureDependentReactions
(
    gasThermoPhysics,
    ArrheniusReactionRate,
    TroeFallOffFunction
)

makePressureDependentReactions
(
    gasThermoPhysics,
    ArrheniusReactionRate,
    SRIFallOffFunction
)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
