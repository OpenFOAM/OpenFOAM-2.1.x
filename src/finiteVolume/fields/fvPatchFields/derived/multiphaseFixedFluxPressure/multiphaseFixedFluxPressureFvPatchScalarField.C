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

#include "multiphaseFixedFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseFixedFluxPressureFvPatchScalarField::
multiphaseFixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    phiHbyAName_("phiHbyA"),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::multiphaseFixedFluxPressureFvPatchScalarField::
multiphaseFixedFluxPressureFvPatchScalarField
(
    const multiphaseFixedFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    phiHbyAName_(ptf.phiHbyAName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::multiphaseFixedFluxPressureFvPatchScalarField::
multiphaseFixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    phiHbyAName_(dict.lookupOrDefault<word>("phiHbyA", "phiHbyA")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::multiphaseFixedFluxPressureFvPatchScalarField::
multiphaseFixedFluxPressureFvPatchScalarField
(
    const multiphaseFixedFluxPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    phiHbyAName_(wbppsf.phiHbyAName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_)
{}


Foam::multiphaseFixedFluxPressureFvPatchScalarField::
multiphaseFixedFluxPressureFvPatchScalarField
(
    const multiphaseFixedFluxPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    phiHbyAName_(wbppsf.phiHbyAName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseFixedFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phiHbyA =
        db().lookupObject<surfaceScalarField>(phiHbyAName_);

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    fvsPatchField<scalar> phiHbyAp =
        patch().patchField<surfaceScalarField, scalar>(phiHbyA);

    fvsPatchField<scalar> phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        phip /= rhop;
    }

    const fvsPatchField<scalar>& Dpp =
        patch().lookupPatchField<surfaceScalarField, scalar>("Dp");

    gradient() = (phiHbyAp - phip)/patch().magSf()/Dpp;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::multiphaseFixedFluxPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phiHbyA", "phiHbyA", phiHbyAName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    gradient().writeEntry("gradient", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        multiphaseFixedFluxPressureFvPatchScalarField
    );
}

// ************************************************************************* //
