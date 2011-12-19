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

#include "waveSurfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    zetaName_("zeta"),
    zeta0_(p.size(), vector::zero),
    curTimeIndex_(-1)
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    zetaName_(dict.lookupOrDefault<word>("zeta", "zeta")),
    zeta0_(p.size(), vector::zero),
    curTimeIndex_(-1)
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    zetaName_(ptf.zetaName_),
    zeta0_(ptf.zeta0_),
    curTimeIndex_(-1)
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& wspsf
)
:
    fixedValueFvPatchScalarField(wspsf),
    phiName_(wspsf.phiName_),
    rhoName_(wspsf.rhoName_),
    zetaName_(wspsf.zetaName_),
    zeta0_(wspsf.zeta0_),
    curTimeIndex_(wspsf.curTimeIndex_)
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& wspsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wspsf, iF),
    phiName_(wspsf.phiName_),
    rhoName_(wspsf.rhoName_),
    zetaName_(wspsf.zetaName_),
    zeta0_(wspsf.zeta0_),
    curTimeIndex_(wspsf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveSurfacePressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    zeta0_.autoMap(m);
}


void Foam::waveSurfacePressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const waveSurfacePressureFvPatchScalarField& wspsf =
        refCast<const waveSurfacePressureFvPatchScalarField>(ptf);

    zeta0_.rmap(wspsf.zeta0_, addr);
}


void Foam::waveSurfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar dt = db().time().deltaTValue();
    const scalar timeI = db().time().timeIndex();
    const scalar patchI = patch().index();

    volVectorField& zeta =
        const_cast<volVectorField&>
        (
            db().lookupObject<volVectorField>(zetaName_)
        );

    vectorField& zetap = zeta.boundaryField()[patchI];

    if (curTimeIndex_ != timeI)
    {
        zeta0_ = zetap;
        curTimeIndex_ = timeI;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const scalarField& phip = phi.boundaryField()[patchI];

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    tmp<vectorField> nf(patch().nf());

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        zetap = zeta0_ + nf()*dt*phip/patch().magSf();

        operator==(-g.value() & zetap);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const scalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        zetap = zeta0_ + nf()*dt*phip/rhop/patch().magSf();

        operator==(-rhop*(g.value() & zetap));
    }
    else
    {
        FatalErrorIn
        (
            "waveSurfacePressureFvPatchScalarField::updateCoeffs()"
        )
            << "dimensions of phi are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    Info<< "min/max mag(zetap) = " << min(zetap & nf()) << ", "
        << max(zetap & nf()) << endl;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::waveSurfacePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "zeta", "zeta", zetaName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waveSurfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
