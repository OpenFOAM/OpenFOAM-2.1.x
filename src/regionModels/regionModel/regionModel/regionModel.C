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

#include "regionModel.H"
#include "fvMesh.H"
#include "Time.H"
#include "mappedWallPolyPatch.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionModel, 0);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::regionModel::constructMeshObjects()
{
    // construct region mesh
    if (!time_.foundObject<fvMesh>(regionName_))
    {
        regionMeshPtr_.reset
        (
            new fvMesh
            (
                IOobject
                (
                    regionName_,
                    time_.timeName(),
                    time_,
                    IOobject::MUST_READ
                )
            )
        );
    }
}


void Foam::regionModels::regionModel::constructMeshObjects
(
    const dictionary& dict
)
{
    // construct region mesh
    if (!time_.foundObject<fvMesh>(regionName_))
    {
        regionMeshPtr_.reset
        (
            new fvMesh
            (
                IOobject
                (
                    regionName_,
                    time_.timeName(),
                    time_,
                    IOobject::MUST_READ
                )
            )
        );
    }
}


void Foam::regionModels::regionModel::initialise()
{
    if (debug)
    {
        Pout<< "regionModel::initialise()" << endl;
    }

    label nBoundaryFaces = 0;
    DynamicList<label> primaryPatchIDs;
    DynamicList<label> intCoupledPatchIDs;

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    forAll(rbm, patchI)
    {
        const polyPatch& regionPatch = rbm[patchI];
        if (isA<mappedWallPolyPatch>(regionPatch))
        {
            if (debug)
            {
                Pout<< "found " << mappedWallPolyPatch::typeName
                    <<  " " << regionPatch.name() << endl;
            }

            intCoupledPatchIDs.append(patchI);

            nBoundaryFaces += regionPatch.faceCells().size();

            const mappedPatchBase& mapPatch =
                refCast<const mappedPatchBase>(regionPatch);

            if (primaryMesh_.time().foundObject<polyMesh>(mapPatch.sampleRegion()))
            {
                const label primaryPatchI = mapPatch.samplePolyPatch().index();
                primaryPatchIDs.append(primaryPatchI);
            }
        }
    }

    primaryPatchIDs_.transfer(primaryPatchIDs);
    intCoupledPatchIDs_.transfer(intCoupledPatchIDs);

    if (nBoundaryFaces == 0)
    {
        WarningIn("regionModel::initialise()")
            << "Region model has no mapped boundary conditions - transfer "
            << "between regions will not be possible" << endl;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regionModels::regionModel::read()
{
    if (regIOobject::read())
    {
        if (active_)
        {
            if (const dictionary* dictPtr = subDictPtr(modelName_ + "Coeffs"))
            {
                coeffs_ <<= *dictPtr;
            }

            infoOutput_.readIfPresent("infoOutput", *this);
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::regionModels::regionModel::read(const dictionary& dict)
{
    if (active_)
    {
        if (const dictionary* dictPtr = dict.subDictPtr(modelName_ + "Coeffs"))
        {
            coeffs_ <<= *dictPtr;
        }

        infoOutput_.readIfPresent("infoOutput", dict);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionModel::regionModel(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "regionModelProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(false),
    infoOutput_(false),
    modelName_("none"),
    regionMeshPtr_(NULL),
    coeffs_(dictionary::null),
    primaryPatchIDs_(),
    intCoupledPatchIDs_(),
    regionName_("none")
{}


Foam::regionModels::regionModel::regionModel
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    IOdictionary
    (
        IOobject
        (
            regionType + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(lookup("active")),
    infoOutput_(true),
    modelName_(modelName),
    regionMeshPtr_(NULL),
    coeffs_(subOrEmptyDict(modelName + "Coeffs")),
    primaryPatchIDs_(),
    intCoupledPatchIDs_(),
    regionName_(lookup("regionName"))
{
    if (active_)
    {
        constructMeshObjects();
        initialise();

        if (readFields)
        {
            read();
        }
    }
}


Foam::regionModels::regionModel::regionModel
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    const dictionary& dict,
    bool readFields
)
:
    IOdictionary
    (
        IOobject
        (
            regionType,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        dict
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(dict.lookup("active")),
    infoOutput_(false),
    modelName_(modelName),
    regionMeshPtr_(NULL),
    coeffs_(dict.subDict(modelName + "Coeffs")),
    primaryPatchIDs_(),
    intCoupledPatchIDs_(),
    regionName_(dict.lookup("regionName"))
{
    if (active_)
    {
        constructMeshObjects(dict);
        initialise();

        if (readFields)
        {
            read(dict);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::regionModel::~regionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionModels::regionModel::preEvolveRegion()
{
    // do nothing
}


void Foam::regionModels::regionModel::evolveRegion()
{
    // do nothing
}


void Foam::regionModels::regionModel::evolve()
{
    if (active_)
    {

        Info<< "\nEvolving " << modelName_ << " for region "
            << regionMesh().name() << endl;

        // Update any input information
        //read();

        // Pre-evolve
        preEvolveRegion();

        // Increment the region equations up to the new time level
        evolveRegion();

        // Provide some feedback
        if (infoOutput_)
        {
            Info<< incrIndent;
            info();
            Info<< endl << decrIndent;
        }
    }
}


void Foam::regionModels::regionModel::info() const
{
    // do nothing
}


// ************************************************************************* //
