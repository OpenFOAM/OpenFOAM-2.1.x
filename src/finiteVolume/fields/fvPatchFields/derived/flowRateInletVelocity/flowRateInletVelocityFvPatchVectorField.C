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

#include "flowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_().clone().ptr()),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(DataEntry<scalar>::New("flowRate", dict)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (dict.found("value") || (phiName_ != "none"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else 
    {
        evaluate(Pstream::blocking);
    }
}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_().clone().ptr()),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::flowRateInletVelocityFvPatchVectorField::
flowRateInletVelocityFvPatchVectorField
(
    const flowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_().clone().ptr()),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();

    // a simpler way of doing this would be nice
    const scalar avgU = -flowRate_->value(t)/gSum(patch().magSf());

    tmp<vectorField> n = patch().nf();


    if (phiName_ == "none")
    {
        // volumetric flow-rate
        operator==(n*avgU);
    }
    else
    {
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            // volumetric flow-rate
            operator==(n*avgU);
        }
        else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            if (rhoName_ == "none")
            {
                // volumetric flow-rate
                operator==(n*avgU);
            }
            else
            {
                // mass flow-rate
                const fvPatchField<scalar>& rhop =
                    patch().lookupPatchField<volScalarField, scalar>(rhoName_);

                operator==(n*avgU/rhop);
            }
        }
        else
        {
            FatalErrorIn
            (
                "flowRateInletVelocityFvPatchVectorField::updateCoeffs()"
            )   << "dimensions of " << phiName_ << " are incorrect" << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << nl << exit(FatalError);
        }
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::flowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       flowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
