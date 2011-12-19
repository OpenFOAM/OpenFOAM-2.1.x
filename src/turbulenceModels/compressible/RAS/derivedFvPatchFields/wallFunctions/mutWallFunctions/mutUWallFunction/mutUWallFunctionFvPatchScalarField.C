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

#include "mutUWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mutUWallFunctionFvPatchScalarField::calcYPlus
(
    const scalarField& magUp
) const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patchI];
    const fvPatchScalarField& rhow = rasModel.rho().boundaryField()[patchI];
    const fvPatchScalarField& muw = rasModel.mu().boundaryField()[patchI];

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus();

    forAll(yPlus, faceI)
    {
        scalar kappaRe = kappa_*magUp[faceI]*y[faceI]/(muw[faceI]/rhow[faceI]);

        scalar yp = yPlusLam_;
        scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E_*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10);

        yPlus[faceI] = max(0.0, yp);
    }

    return tyPlus;
}


tmp<scalarField> mutUWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const fvPatchScalarField& muw = rasModel.mu().boundaryField()[patchI];

    tmp<scalarField> tyPlus = calcYPlus(magUp);
    scalarField& yPlus = tyPlus();

    tmp<scalarField> tmutw(new scalarField(patch().size(), 0.0));
    scalarField& mutw = tmutw();

    forAll(yPlus, faceI)
    {
        if (yPlus[faceI] > yPlusLam_)
        {
            mutw[faceI] =
                muw[faceI]*(yPlus[faceI]*kappa_/log(E_*yPlus[faceI]) - 1.0);
        }
    }

    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(p, iF)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const mutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutkWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const mutUWallFunctionFvPatchScalarField& sawfpsf
)
:
    mutkWallFunctionFvPatchScalarField(sawfpsf)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const mutUWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(sawfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> mutUWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    return calcYPlus(magUp);
}


void mutUWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutUWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
