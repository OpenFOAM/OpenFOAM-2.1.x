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

#include "htcConvFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

htcConvFvPatchScalarField::htcConvFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    L_(1.0)
{}


htcConvFvPatchScalarField::htcConvFvPatchScalarField
(
    const htcConvFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    L_(ptf.L_)
{}


htcConvFvPatchScalarField::htcConvFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    L_(readScalar(dict.lookup("L")))
{}


htcConvFvPatchScalarField::htcConvFvPatchScalarField
(
    const htcConvFvPatchScalarField& htcpsf
)
:
    fixedValueFvPatchScalarField(htcpsf),
    L_(htcpsf.L_)
{}


htcConvFvPatchScalarField::htcConvFvPatchScalarField
(
    const htcConvFvPatchScalarField& htcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(htcpsf, iF),
    L_(htcpsf.L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void htcConvFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField alphaEffw = rasModel.alphaEff()().boundaryField()[patchI];
    const scalarField& muw = rasModel.mu().boundaryField()[patchI];
    const scalarField& rhow = rasModel.rho().boundaryField()[patchI];
    const vectorField& Uc = rasModel.U();
    const vectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField& Tw = rasModel.thermo().T().boundaryField()[patchI];
    const scalarField Cpw(rasModel.thermo().Cp(Tw, patchI));

    const scalarField kappaw(Cpw*alphaEffw);
    const scalarField Pr(muw*Cpw/kappaw);

    scalarField& htc = *this;
    forAll(htc, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar Re = rhow[faceI]*mag(Uc[faceCellI] - Uw[faceI])*L_/muw[faceI];

        if (Re < 5.0E+05)
        {
            htc[faceI] = 0.664*sqrt(Re)*cbrt(Pr[faceI])*kappaw[faceI]/L_;
        }
        else
        {
            htc[faceI] = 0.037*pow(Re, 0.8)*cbrt(Pr[faceI])*kappaw[faceI]/L_;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void htcConvFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("L") << L_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    htcConvFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
