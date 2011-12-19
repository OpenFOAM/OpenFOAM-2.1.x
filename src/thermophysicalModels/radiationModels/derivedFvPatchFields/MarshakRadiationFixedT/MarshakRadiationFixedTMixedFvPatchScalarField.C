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

#include "MarshakRadiationFixedTMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvc.H"
#include "radiationModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    Trad_(p.size())
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const MarshakRadiationFixedTMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    radiationCoupledBase
    (
        p,
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    Trad_(ptf.Trad_, mapper)
{}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, dict),
    Trad_("Trad", dict, p.size())
{
    // refValue updated on each call to updateCoeffs()
    refValue() = 4.0*constant::physicoChemical::sigma.value()*pow4(Trad_);

    // zero gradient
    refGrad() = 0.0;

    valueFraction() = 1.0;

    fvPatchScalarField::operator=(refValue());
}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const MarshakRadiationFixedTMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    Trad_(ptf.Trad_)
{}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const MarshakRadiationFixedTMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    radiationCoupledBase
    (
        ptf.patch(),
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    Trad_(ptf.Trad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    Trad_.autoMap(m);
}


void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const MarshakRadiationFixedTMixedFvPatchScalarField& mrptf =
        refCast<const MarshakRadiationFixedTMixedFvPatchScalarField>(ptf);

    Trad_.rmap(mrptf.Trad_, addr);
}


void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Re-calc reference value
    refValue() = 4.0*constant::physicoChemical::sigma.value()*pow4(Trad_);

    // Diffusion coefficient - created by radiation model's ::updateCoeffs()
    const scalarField& gamma =
        patch().lookupPatchField<volScalarField, scalar>("gammaRad");

    const scalarField temissivity = emissivity();

    const scalarField Ep(temissivity/(2.0*(scalar(2.0) - temissivity)));

    // Set value fraction
    valueFraction() = 1.0/(1.0 + gamma*patch().deltaCoeffs()/Ep);

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    radiationCoupledBase::write(os);
    Trad_.writeEntry("Trad", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        MarshakRadiationFixedTMixedFvPatchScalarField
    );
}


// ************************************************************************* //
