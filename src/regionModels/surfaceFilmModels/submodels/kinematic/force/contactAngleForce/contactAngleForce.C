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

#include "contactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "unitConversion.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactAngleForce, 0);
addToRunTimeSelectionTable(force, contactAngleForce, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactAngleForce::contactAngleForce
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    force(typeName, owner, dict),
    deltaWet_(readScalar(coeffs_.lookup("deltaWet"))),
    Ccf_(readScalar(coeffs_.lookup("Ccf"))),
    rndGen_(label(0), -1),
    distribution_
    (
        distributionModels::distributionModel::New
        (
            coeffs_.subDict("contactAngleDistribution"),
            rndGen_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

contactAngleForce::~contactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> contactAngleForce::correct(volVectorField& U)
{
    tmp<volVectorField> tForce
    (
        new volVectorField
        (
            IOobject
            (
                "contactForce",
                owner_.time().timeName(),
                owner_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner_.regionMesh(),
            dimensionedVector("zero", dimForce/dimArea, vector::zero)
        )
    );

    vectorField& force = tForce().internalField();

    const labelUList& own = owner_.regionMesh().owner();
    const labelUList& nbr = owner_.regionMesh().neighbour();

    const scalarField& magSf = owner_.magSf();

    const volScalarField& delta = owner_.delta();
    const volScalarField& sigma = owner_.sigma();

    volScalarField alpha
    (
        "alpha",
        pos(delta - dimensionedScalar("deltaWet", dimLength, deltaWet_))
    );
    volVectorField gradAlpha(fvc::grad(alpha));


    scalarField nHits(owner_.regionMesh().nCells(), 0.0);

    forAll(nbr, faceI)
    {
        const label cellO = own[faceI];
        const label cellN = nbr[faceI];

        label cellI = -1;
        if ((delta[cellO] > deltaWet_) && (delta[cellN] < deltaWet_))
        {
            cellI = cellO;
        }
        else if ((delta[cellO] < deltaWet_) && (delta[cellN] > deltaWet_))
        {
            cellI = cellN;
        }

        if (cellI != -1)
        {
//            const scalar dx = Foam::sqrt(magSf[cellI]);
            // bit of a cheat, but ok for regular meshes
            const scalar dx = owner_.regionMesh().deltaCoeffs()[faceI];
            const vector n =
                gradAlpha[cellI]/(mag(gradAlpha[cellI]) + ROOTVSMALL);
            scalar theta = cos(degToRad(distribution_->sample()));
            force[cellI] += Ccf_*n*sigma[cellI]*(1.0 - theta)/dx;
            nHits[cellI]++;
        }
    }

    forAll(delta.boundaryField(), patchI)
    {
        const fvPatchField<scalar>& df = delta.boundaryField()[patchI];
        const scalarField& dx = df.patch().deltaCoeffs();
        const labelUList& faceCells = df.patch().faceCells();

        forAll(df, faceI)
        {
            label cellO = faceCells[faceI];

            if ((delta[cellO] > deltaWet_) && (df[faceI] < deltaWet_))
            {
                const vector n =
                    gradAlpha[cellO]/(mag(gradAlpha[cellO]) + ROOTVSMALL);
                scalar theta = cos(degToRad(distribution_->sample()));
                force[cellO] += Ccf_*n*sigma[cellO]*(1.0 - theta)/dx[faceI];
                nHits[cellO]++;
            }
        }
    }

    force /= (max(nHits, scalar(1.0))*magSf);
    tForce().correctBoundaryConditions();

    if (owner_.regionMesh().time().outputTime())
    {
        tForce().write();
    }

    tmp<fvVectorMatrix>
        tfvm(new fvVectorMatrix(U, dimForce/dimArea*dimVolume));

    tfvm() += tForce;

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
