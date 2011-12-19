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

#include "filmPyrolysisVelocityCoupledFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "pyrolysisModel.H"
#include "surfaceFilmModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    deltaWet_(1e-6)
{}


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const filmPyrolysisVelocityCoupledFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    deltaWet_(ptf.deltaWet_)
{}


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    deltaWet_(dict.lookupOrDefault<scalar>("deltaWet", 1e-6))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const filmPyrolysisVelocityCoupledFvPatchVectorField& fpvpvf
)
:
    fixedValueFvPatchVectorField(fpvpvf),
    phiName_(fpvpvf.phiName_),
    rhoName_(fpvpvf.rhoName_),
    deltaWet_(fpvpvf.deltaWet_)
{}


Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::
filmPyrolysisVelocityCoupledFvPatchVectorField
(
    const filmPyrolysisVelocityCoupledFvPatchVectorField& fpvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fpvpvf, iF),
    phiName_(fpvpvf.phiName_),
    rhoName_(fpvpvf.rhoName_),
    deltaWet_(fpvpvf.deltaWet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    typedef regionModels::surfaceFilmModels::surfaceFilmModel filmModelType;
    typedef regionModels::pyrolysisModels::pyrolysisModel pyrModelType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    bool filmOk =
        db().objectRegistry::foundObject<filmModelType>
        (
            "surfaceFilmProperties"
        );


    bool pyrOk =
        db().objectRegistry::foundObject<pyrModelType>
        (
            "pyrolysisProperties"
        );

    if (!filmOk || !pyrOk)
    {
        // do nothing on construction - film model doesn't exist yet
        return;
    }

    vectorField& Up = *this;

    const label patchI = patch().index();

    // Retrieve film model
    const filmModelType& filmModel =
        db().objectRegistry::lookupObject<filmModelType>
        (
            "surfaceFilmProperties"
        );

    const label filmPatchI = filmModel.regionPatchID(patchI);

    scalarField deltaFilm = filmModel.delta().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, deltaFilm);

    vectorField UFilm = filmModel.Us().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, UFilm);


    // Retrieve pyrolysis model
    const pyrModelType& pyrModel =
        db().objectRegistry::lookupObject<pyrModelType>
        (
            "pyrolysisProperties"
        );

    const label pyrPatchI = pyrModel.regionPatchID(patchI);

    scalarField phiPyr = pyrModel.phiGas().boundaryField()[pyrPatchI];
    pyrModel.toPrimary(pyrPatchI, phiPyr);


    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // do nothing
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        phiPyr /= rhop;
    }
    else
    {
        FatalErrorIn
        (
            "filmPyrolysisVelocityCoupledFvPatchVectorField::updateCoeffs()"
        )   << "Unable to process flux field phi with dimensions "
            << phi.dimensions() << nl
            << "    on patch " << patch().name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    const scalarField UAvePyr(-phiPyr/patch().magSf());
    const vectorField& nf = patch().nf();

    forAll(deltaFilm, i)
    {
        if (deltaFilm[i] > deltaWet_)
        {
            // velocity set by film
            Up[i] = UFilm[i];
        }
        else
        {
            // velocity set by pyrolysis model
            Up[i] = UAvePyr[i]*nf[i];
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::filmPyrolysisVelocityCoupledFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    os.writeKeyword("deltaWet") << deltaWet_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        filmPyrolysisVelocityCoupledFvPatchVectorField
    );
}


// ************************************************************************* //
