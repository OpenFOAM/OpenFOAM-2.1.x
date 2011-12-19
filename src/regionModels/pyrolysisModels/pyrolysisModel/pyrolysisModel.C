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

#include "pyrolysisModel.H"
#include "fvMesh.H"
#include "mappedFieldFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pyrolysisModel, 0);
defineRunTimeSelectionTable(pyrolysisModel, mesh);
defineRunTimeSelectionTable(pyrolysisModel, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void pyrolysisModel::constructMeshObjects()
{
    // construct filmDelta field if coupled to film model
    if (filmCoupled_)
    {
        filmDeltaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "filmDelta",
                    time_.timeName(),
                    regionMesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                regionMesh()
            )
        );

        const volScalarField& filmDelta = filmDeltaPtr_();

        bool foundCoupledPatch = false;
        forAll(filmDelta.boundaryField(), patchI)
        {
            const fvPatchField<scalar>& fvp = filmDelta.boundaryField()[patchI];
            if (isA<mappedFieldFvPatchField<scalar> >(fvp))
            {
                foundCoupledPatch = true;
                break;
            }
        }

        if (!foundCoupledPatch)
        {
            WarningIn("void pyrolysisModels::constructMeshObjects()")
                << "filmCoupled flag set to true, but no "
                << mappedFieldFvPatchField<scalar>::typeName
                << " patches found on " << filmDelta.name() << " field"
                << endl;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void pyrolysisModel::readPyrolysisControls()
{
    filmCoupled_ = readBool(coeffs_.lookup("filmCoupled"));
    reactionDeltaMin_ =
        coeffs_.lookupOrDefault<scalar>("reactionDeltaMin", 0.0);
}


bool pyrolysisModel::read()
{
    if (regionModel1D::read())
    {
        readPyrolysisControls();
        return true;
    }
    else
    {
        return false;
    }
}


bool pyrolysisModel::read(const dictionary& dict)
{
    if (regionModel1D::read(dict))
    {
        readPyrolysisControls();
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pyrolysisModel::pyrolysisModel(const fvMesh& mesh)
:
    regionModel1D(mesh),
    filmCoupled_(false),
    filmDeltaPtr_(NULL),
    reactionDeltaMin_(0.0)
{}


pyrolysisModel::pyrolysisModel(const word& modelType, const fvMesh& mesh)
:
    regionModel1D(mesh, "pyrolysis", modelType),
    filmCoupled_(false),
    filmDeltaPtr_(NULL),
    reactionDeltaMin_(0.0)
{
    if (active_)
    {
        read();
        constructMeshObjects();
    }
}


pyrolysisModel::pyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regionModel1D(mesh, "pyrolysis", modelType, dict),
    filmCoupled_(false),
    filmDeltaPtr_(NULL),
    reactionDeltaMin_(0.0)
{
    if (active_)
    {
        read(dict);
        constructMeshObjects();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pyrolysisModel::~pyrolysisModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar pyrolysisModel::addMassSources
(
    const label patchI,
    const label faceI
)
{
    return 0.0;
}


void pyrolysisModel::preEvolveRegion()
{
    if (filmCoupled_)
    {
        filmDeltaPtr_->correctBoundaryConditions();
    }
}


scalar pyrolysisModel::solidRegionDiffNo() const
{
    return VSMALL;
}


scalar pyrolysisModel::maxDiff() const
{
    return GREAT;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
