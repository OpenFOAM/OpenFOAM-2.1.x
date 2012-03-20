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

#include "faceSource.H"
#include "fvMesh.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "sampledSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fieldValues::faceSource, 0);

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::fieldValues::faceSource::sourceType,
        3
    >::names[] =
    {
        "faceZone",
        "patch",
        "sampledSurface"
    };


    template<>
    const char* Foam::NamedEnum
    <
        Foam::fieldValues::faceSource::operationType,
        11
    >::names[] =
    {
        "none",
        "sum",
        "average",
        "weightedAverage",
        "areaAverage",
        "areaIntegrate",
        "min",
        "max",
        "CoV",
        "areaNormalAverage",
        "areaNormalIntegrate"
    };

}


const Foam::NamedEnum<Foam::fieldValues::faceSource::sourceType, 3>
    Foam::fieldValues::faceSource::sourceTypeNames_;

const Foam::NamedEnum<Foam::fieldValues::faceSource::operationType, 11>
    Foam::fieldValues::faceSource::operationTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::faceSource::setFaceZoneFaces()
{
    label zoneId = mesh().faceZones().findZoneID(sourceName_);

    if (zoneId < 0)
    {
        FatalErrorIn("faceSource::faceSource::setFaceZoneFaces()")
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Unknown face zone name: " << sourceName_
            << ". Valid face zones are: " << mesh().faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh().faceZones()[zoneId];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        label faceI = fZone[i];

        label faceId = -1;
        label facePatchId = -1;
        if (mesh().isInternalFace(faceI))
        {
            faceId = faceI;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh().boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh().boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(faceI);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = faceI - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }

    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    if (debug)
    {
        Pout<< "Original face zone size = " << fZone.size()
            << ", new size = " << count << endl;
    }
}


void Foam::fieldValues::faceSource::setPatchFaces()
{
    const label patchId = mesh().boundaryMesh().findPatchID(sourceName_);

    if (patchId < 0)
    {
        FatalErrorIn("faceSource::constructFaceAddressing()")
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Unknown patch name: " << sourceName_
            << ". Valid patch names are: "
            << mesh().boundaryMesh().names() << nl
            << exit(FatalError);
    }

    const polyPatch& pp = mesh().boundaryMesh()[patchId];

    label nFaces = pp.size();
    if (isA<cyclicPolyPatch>(pp))
    {
        nFaces /= 2;
    }
    else if (isA<emptyPolyPatch>(pp))
    {
        nFaces = 0;
    }

    faceId_.setSize(nFaces);
    facePatchId_.setSize(nFaces);
    faceSign_.setSize(nFaces);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    forAll(faceId_, faceI)
    {
        faceId_[faceI] = faceI;
        facePatchId_[faceI] = patchId;
        faceSign_[faceI] = 1;
    }
}


void Foam::fieldValues::faceSource::sampledSurfaceFaces(const dictionary& dict)
{
    surfacePtr_ = sampledSurface::New
    (
        name_,
        mesh(),
        dict.subDict("sampledSurfaceDict")
    );
    surfacePtr_().update();
    nFaces_ = returnReduce(surfacePtr_().faces().size(), sumOp<label>());
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::faceSource::initialise(const dictionary& dict)
{
    switch (source_)
    {
        case stFaceZone:
        {
            setFaceZoneFaces();
            break;
        }
        case stPatch:
        {
            setPatchFaces();
            break;
        }
        case stSampledSurface:
        {
            sampledSurfaceFaces(dict);
            break;
        }
        default:
        {
            FatalErrorIn("faceSource::initialise()")
                << type() << " " << name_ << ": "
                << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                << nl << "    Unknown source type. Valid source types are:"
                << sourceTypeNames_.sortedToc() << nl << exit(FatalError);
        }
    }

    if (nFaces_ == 0)
    {
        WarningIn
        (
            "Foam::fieldValues::faceSource::initialise(const dictionary&)"
        )
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Source has no faces - deactivating" << endl;

        active_ = false;
        return;
    }

    scalar totalArea;

    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
        totalArea = gSum(surfacePtr_().magSf());
    }
    else
    {
        totalArea = gSum(filterField(mesh().magSf(), false));
    }

    Info<< type() << " " << name_ << ":" << nl
        << "    total faces  = " << nFaces_
        << nl
        << "    total area   = " << totalArea
        << nl;

    if (operation_ == opWeightedAverage)
    {
        dict.lookup("weightField") >> weightFieldName_;
        Info<< "    weight field = " << weightFieldName_;
    }

    Info<< nl << endl;
}


void Foam::fieldValues::faceSource::writeFileHeader()
{
    if (outputFilePtr_.valid())
    {
        outputFilePtr_()
            << "# Source : " << sourceTypeNames_[source_] << " "
            << sourceName_ <<  nl << "# Faces  : " << nFaces_ << nl
            << "# Time" << tab << "sum(magSf)";

        forAll(fields_, i)
        {
            outputFilePtr_()
                << tab << operationTypeNames_[operation_]
                << "(" << fields_[i] << ")";
        }

        outputFilePtr_() << endl;
    }
}


template<>
Foam::vector Foam::fieldValues::faceSource::processValues
(
    const Field<vector>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opAreaNormalAverage:
        {
            scalar result = sum(values&Sf)/sum(mag(Sf));
            return vector(result, 0.0, 0.0);
        }
        case opAreaNormalIntegrate:
        {
            scalar result = sum(values&Sf);
            return vector(result, 0.0, 0.0);
        }
        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValues::faceSource::faceSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    fieldValue(name, obr, dict, loadFromFiles),
    source_(sourceTypeNames_.read(dict.lookup("source"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldName_("none"),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceSign_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValues::faceSource::~faceSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValues::faceSource::read(const dictionary& dict)
{
    fieldValue::read(dict);

    if (active_)
    {
        initialise(dict);
    }
}


void Foam::fieldValues::faceSource::write()
{
    fieldValue::write();

    if (active_)
    {
        scalar totalArea;

        if (surfacePtr_.valid())
        {
            surfacePtr_().update();
            totalArea = gSum(surfacePtr_().magSf());
        }
        else
        {
            totalArea = gSum(filterField(mesh().magSf(), false));
        }


        if (Pstream::master())
        {
            outputFilePtr_() << obr_.time().value() << tab << totalArea;
        }

        forAll(fields_, i)
        {
            writeValues<scalar>(fields_[i]);
            writeValues<vector>(fields_[i]);
            writeValues<sphericalTensor>(fields_[i]);
            writeValues<symmTensor>(fields_[i]);
            writeValues<tensor>(fields_[i]);
        }

        if (Pstream::master())
        {
            outputFilePtr_()<< endl;
        }

        if (log_)
        {
            Info<< endl;
        }
    }
}


// ************************************************************************* //
