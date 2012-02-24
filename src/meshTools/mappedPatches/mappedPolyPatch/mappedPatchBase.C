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

#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearch.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Random.H"
#include "treeDataFace.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "Time.H"
#include "mapDistribute.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBase, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::mappedPatchBase::sampleMode,
        4
    >::names[] =
    {
        "nearestCell",
        "nearestPatchFace",
        "nearestPatchFaceAMI",
        "nearestFace"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::mappedPatchBase::offsetMode,
        3
    >::names[] =
    {
        "uniform",
        "nonuniform",
        "normal"
    };
}


const Foam::NamedEnum<Foam::mappedPatchBase::sampleMode, 4>
    Foam::mappedPatchBase::sampleModeNames_;

const Foam::NamedEnum<Foam::mappedPatchBase::offsetMode, 3>
    Foam::mappedPatchBase::offsetModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mappedPatchBase::collectSamples
(
    pointField& samples,
    labelList& patchFaceProcs,
    labelList& patchFaces,
    pointField& patchFc
) const
{

    // Collect all sample points and the faces they come from.
    List<pointField> globalFc(Pstream::nProcs());
    List<pointField> globalSamples(Pstream::nProcs());
    labelListList globalFaces(Pstream::nProcs());

    globalFc[Pstream::myProcNo()] = patch_.faceCentres();
    globalSamples[Pstream::myProcNo()] = samplePoints();
    globalFaces[Pstream::myProcNo()] = identity(patch_.size());

    // Distribute to all processors
    Pstream::gatherList(globalSamples);
    Pstream::scatterList(globalSamples);
    Pstream::gatherList(globalFaces);
    Pstream::scatterList(globalFaces);
    Pstream::gatherList(globalFc);
    Pstream::scatterList(globalFc);

    // Rework into straight list
    samples = ListListOps::combine<pointField>
    (
        globalSamples,
        accessOp<pointField>()
    );
    patchFaces = ListListOps::combine<labelList>
    (
        globalFaces,
        accessOp<labelList>()
    );
    patchFc = ListListOps::combine<pointField>
    (
        globalFc,
        accessOp<pointField>()
    );

    patchFaceProcs.setSize(patchFaces.size());
    labelList nPerProc
    (
        ListListOps::subSizes
        (
            globalFaces,
            accessOp<labelList>()
        )
    );
    label sampleI = 0;
    forAll(nPerProc, procI)
    {
        for (label i = 0; i < nPerProc[procI]; i++)
        {
            patchFaceProcs[sampleI++] = procI;
        }
    }
}


// Find the processor/cell containing the samples. Does not account
// for samples being found in two processors.
void Foam::mappedPatchBase::findSamples
(
    const pointField& samples,
    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    // Lookup the correct region
    const polyMesh& mesh = sampleMesh();

    // All the info for nearest. Construct to miss
    List<nearInfo> nearest(samples.size());

    switch (mode_)
    {
        case NEARESTCELL:
        {
            if (samplePatch_.size() && samplePatch_ != "none")
            {
                FatalErrorIn
                (
                    "mappedPatchBase::findSamples(const pointField&,"
                    " labelList&, labelList&, pointField&) const"
                )   << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            // Octree based search engine. Uses min tetDecomp so force
            // calculation
            (void)mesh.tetBasePtIs();
            meshSearch meshSearchEngine(mesh);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label cellI = meshSearchEngine.findCell(sample);

                if (cellI == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& cc = mesh.cellCentres()[cellI];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        cc,
                        cellI
                    );
                    nearest[sampleI].second().first() = magSqr(cc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            Random rndGen(123456);

            const polyPatch& pp = samplePolyPatch();

            if (pp.empty())
            {
                forAll(samples, sampleI)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            else
            {
                // patch faces
                const labelList patchFaces(identity(pp.size()) + pp.start());

                treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend
                    (
                        rndGen,
                        1E-4
                    )
                );
                patchBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
                patchBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

                indexedOctree<treeDataFace> boundaryTree
                (
                    treeDataFace    // all information needed to search faces
                    (
                        false,      // do not cache bb
                        mesh,
                        patchFaces  // boundary faces only
                    ),
                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    pointIndexHit& nearInfo = nearest[sampleI].first();
                    nearInfo = boundaryTree.findNearest
                    (
                        sample,
                        magSqr(patchBb.span())
                    );

                    if (!nearInfo.hit())
                    {
                        nearest[sampleI].second().first() = Foam::sqr(GREAT);
                        nearest[sampleI].second().second() =
                            Pstream::myProcNo();
                    }
                    else
                    {
                        point fc(pp[nearInfo.index()].centre(pp.points()));
                        nearInfo.setPoint(fc);
                        nearest[sampleI].second().first() = magSqr(fc-sample);
                        nearest[sampleI].second().second() =
                            Pstream::myProcNo();
                    }
                }
            }
            break;
        }

        case NEARESTFACE:
        {
            if (samplePatch_.size() && samplePatch_ != "none")
            {
                FatalErrorIn
                (
                    "mappedPatchBase::findSamples(const pointField&,"
                    " labelList&, labelList&, pointField&) const"
                )   << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            // Octree based search engine
            meshSearch meshSearchEngine(mesh);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label faceI = meshSearchEngine.findNearestFace(sample);

                if (faceI == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& fc = mesh.faceCentres()[faceI];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        fc,
                        faceI
                    );
                    nearest[sampleI].second().first() = magSqr(fc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        case NEARESTPATCHFACEAMI:
        {
            // nothing to do here
            return;
        }

        default:
        {
            FatalErrorIn("mappedPatchBase::findSamples(..)")
                << "problem." << abort(FatalError);
        }
    }


    // Find nearest
    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);

    if (debug)
    {
        Info<< "mappedPatchBase::findSamples on mesh " << sampleRegion_
            << " : " << endl;
        forAll(nearest, sampleI)
        {
            label procI = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< samples[sampleI]
                << " found on processor:" << procI
                << " in local cell/face:" << localI
                << " with cc:" << nearest[sampleI].first().rawPoint() << endl;
        }
    }

    // Check for samples not being found
    forAll(nearest, sampleI)
    {
        if (!nearest[sampleI].first().hit())
        {
            FatalErrorIn
            (
                "mappedPatchBase::findSamples"
                "(const pointField&, labelList&"
                ", labelList&, pointField&)"
            )   << "Did not find sample " << samples[sampleI]
                << " on any processor of region " << sampleRegion_
                << exit(FatalError);
        }
    }


    // Convert back into proc+local index
    sampleProcs.setSize(samples.size());
    sampleIndices.setSize(samples.size());
    sampleLocations.setSize(samples.size());

    forAll(nearest, sampleI)
    {
        sampleProcs[sampleI] = nearest[sampleI].second().second();
        sampleIndices[sampleI] = nearest[sampleI].first().index();
        sampleLocations[sampleI] = nearest[sampleI].first().hitPoint();
    }
}


void Foam::mappedPatchBase::calcMapping() const
{
    if (mapPtr_.valid())
    {
        FatalErrorIn("mappedPatchBase::calcMapping() const")
            << "Mapping already calculated" << exit(FatalError);
    }

    // Do a sanity check
    // Am I sampling my own patch? This only makes sense for a non-zero
    // offset.
    bool sampleMyself =
    (
        mode_ == NEARESTPATCHFACE
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
     && samplePatch_ == patch_.name()
    );

    // Check offset
    vectorField d(samplePoints()-patch_.faceCentres());
    if (sampleMyself && gAverage(mag(d)) <= ROOTVSMALL)
    {
        WarningIn
        (
            "mappedPatchBase::mappedPatchBase\n"
            "(\n"
            "    const polyPatch& pp,\n"
            "    const word& sampleRegion,\n"
            "    const sampleMode mode,\n"
            "    const word& samplePatch,\n"
            "    const vector& offset\n"
            ")\n"
        )   << "Invalid offset " << d << endl
            << "Offset is the vector added to the patch face centres to"
            << " find the patch face supplying the data." << endl
            << "Setting it to " << d
            << " on the same patch, on the same region"
            << " will find the faces themselves which does not make sense"
            << " for anything but testing." << endl
            << "patch_:" << patch_.name() << endl
            << "sampleRegion_:" << sampleRegion_ << endl
            << "mode_:" << sampleModeNames_[mode_] << endl
            << "samplePatch_:" << samplePatch_ << endl
            << "offsetMode_:" << offsetModeNames_[offsetMode_] << endl;
    }


    // Get global list of all samples and the processor and face they come from.
    pointField samples;
    labelList patchFaceProcs;
    labelList patchFaces;
    pointField patchFc;
    collectSamples(samples, patchFaceProcs, patchFaces, patchFc);

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    pointField sampleLocations;
    findSamples(samples, sampleProcs, sampleIndices, sampleLocations);


    // Now we have all the data we need:
    // - where sample originates from (so destination when mapping):
    //   patchFaces, patchFaceProcs.
    // - cell/face sample is in (so source when mapping)
    //   sampleIndices, sampleProcs.

    //forAll(samples, i)
    //{
    //    Info<< i << " need data in region "
    //        << patch_.boundaryMesh().mesh().name()
    //        << " for proc:" << patchFaceProcs[i]
    //        << " face:" << patchFaces[i]
    //        << " at:" << patchFc[i] << endl
    //        << "Found data in region " << sampleRegion_
    //        << " at proc:" << sampleProcs[i]
    //        << " face:" << sampleIndices[i]
    //        << " at:" << sampleLocations[i]
    //        << nl << endl;
    //}



    if (debug && Pstream::master())
    {
        OFstream str
        (
            patch_.boundaryMesh().mesh().time().path()
          / patch_.name()
          + "_mapped.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cell/faceCentres to file " << str.name() << endl;

        label vertI = 0;

        forAll(patchFc, i)
        {
            meshTools::writeOBJ(str, patchFc[i]);
            vertI++;
            meshTools::writeOBJ(str, sampleLocations[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }


    // Determine schedule.
    mapPtr_.reset(new mapDistribute(sampleProcs, patchFaceProcs));

    // Rework the schedule from indices into samples to cell data to send,
    // face data to receive.

    labelListList& subMap = mapPtr_().subMap();
    labelListList& constructMap = mapPtr_().constructMap();

    forAll(subMap, procI)
    {
        subMap[procI] = UIndirectList<label>
        (
            sampleIndices,
            subMap[procI]
        );
        constructMap[procI] = UIndirectList<label>
        (
            patchFaces,
            constructMap[procI]
        );

        //if (debug)
        //{
        //    Pout<< "To proc:" << procI << " sending values of cells/faces:"
        //        << subMap[procI] << endl;
        //    Pout<< "From proc:" << procI
        //        << " receiving values of patch faces:"
        //        << constructMap[procI] << endl;
        //}
    }

    // Redo constructSize
    mapPtr_().constructSize() = patch_.size();

    if (debug)
    {
        // Check that all elements get a value.
        PackedBoolList used(patch_.size());
        forAll(constructMap, procI)
        {
            const labelList& map = constructMap[procI];

            forAll(map, i)
            {
                label faceI = map[i];

                if (used[faceI] == 0)
                {
                    used[faceI] = 1;
                }
                else
                {
                    FatalErrorIn("mappedPatchBase::calcMapping() const")
                        << "On patch " << patch_.name()
                        << " patchface " << faceI
                        << " is assigned to more than once."
                        << abort(FatalError);
                }
            }
        }
        forAll(used, faceI)
        {
            if (used[faceI] == 0)
            {
                FatalErrorIn("mappedPatchBase::calcMapping() const")
                    << "On patch " << patch_.name()
                    << " patchface " << faceI
                    << " is never assigned to."
                    << abort(FatalError);
            }
        }
    }
}


const Foam::autoPtr<Foam::searchableSurface>& Foam::mappedPatchBase::surfPtr()
const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", patch_.name()));

        const polyMesh& mesh = patch_.boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


void Foam::mappedPatchBase::calcAMI() const
{
    if (AMIPtr_.valid())
    {
        FatalErrorIn("mappedPatchBase::calcAMI() const")
            << "AMI already calculated" << exit(FatalError);
    }

    AMIPtr_.clear();
/*
    const polyPatch& nbr = samplePolyPatch();

//    pointField nbrPoints(samplePoints());
    pointField nbrPoints(nbr.localPoints());

    if (debug)
    {
        OFstream os(patch_.name() + "_neighbourPatch-org.obj");
        meshTools::writeOBJ(os, samplePolyPatch().localFaces(), nbrPoints);
    }

    // transform neighbour patch to local system
    primitivePatch nbrPatch0
    (
        SubList<face>
        (
            nbr.localFaces(),
            nbr.size()
        ),
        nbrPoints
    );

    if (debug)
    {
        OFstream osN(patch_.name() + "_neighbourPatch-trans.obj");
        meshTools::writeOBJ(osN, nbrPatch0, nbrPoints);

        OFstream osO(patch_.name() + "_ownerPatch.obj");
        meshTools::writeOBJ(osO, patch_.localFaces(), patch_.localPoints());
    }
*/
    // Construct/apply AMI interpolation to determine addressing and weights
    AMIPtr_.reset
    (
        new AMIPatchToPatchInterpolation
        (
            patch_,
            samplePolyPatch(), // nbrPatch0,
            surfPtr(),
            faceAreaIntersect::tmMesh,
            AMIReverse_
        )
    );
}


// Hack to read old (List-based) format. See Field.C. The difference
// is only that in case of falling back to old format it expects a non-uniform
// list instead of a single vector.
Foam::tmp<Foam::pointField> Foam::mappedPatchBase::readListOrField
(
    const word& keyword,
    const dictionary& dict,
    const label size
)
{
    tmp<pointField> tfld(new pointField());
    pointField& fld = tfld();

    if (size)
    {
        ITstream& is = dict.lookup(keyword);

        // Read first token
        token firstToken(is);

        if (firstToken.isWord())
        {
            if (firstToken.wordToken() == "uniform")
            {
                fld.setSize(size);
                fld = pTraits<vector>(is);
            }
            else if (firstToken.wordToken() == "nonuniform")
            {
                is >> static_cast<List<vector>&>(fld);
                if (fld.size() != size)
                {
                    FatalIOErrorIn
                    (
                        "mappedPatchBase::readListOrField"
                        "(const word& keyword, const dictionary&, const label)",
                        dict
                    )   << "size " << fld.size()
                        << " is not equal to the given value of " << size
                        << exit(FatalIOError);
                }
            }
            else
            {
                FatalIOErrorIn
                (
                    "mappedPatchBase::readListOrField"
                    "(const word& keyword, const dictionary&, const label)",
                    dict
                )   << "expected keyword 'uniform' or 'nonuniform', found "
                    << firstToken.wordToken()
                    << exit(FatalIOError);
            }
        }
        else
        {
            if (is.version() == 2.0)
            {
                IOWarningIn
                (
                    "mappedPatchBase::readListOrField"
                    "(const word& keyword, const dictionary&, const label)",
                    dict
                )   << "expected keyword 'uniform' or 'nonuniform', "
                       "assuming List format for backwards compatibility."
                       "Foam version 2.0." << endl;

                is.putBack(firstToken);
                is >> static_cast<List<vector>&>(fld);
            }
        }
    }
    return tfld;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp
)
:
    patch_(pp),
    sampleRegion_(patch_.boundaryMesh().mesh().name()),
    mode_(NEARESTPATCHFACE),
    samplePatch_("none"),
    offsetMode_(UNIFORM),
    offset_(vector::zero),
    offsets_(pp.size(), offset_),
    distance_(0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL),
    AMIPtr_(NULL),
    AMIReverse_(false),
    surfPtr_(NULL),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vectorField& offsets
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    offsetMode_(NONUNIFORM),
    offset_(vector::zero),
    offsets_(offsets),
    distance_(0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL),
    AMIPtr_(NULL),
    AMIReverse_(false),
    surfPtr_(NULL),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vector& offset
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    offsetMode_(UNIFORM),
    offset_(offset),
    offsets_(0),
    distance_(0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL),
    AMIPtr_(NULL),
    AMIReverse_(false),
    surfPtr_(NULL),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const scalar distance
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    offsetMode_(NORMAL),
    offset_(vector::zero),
    offsets_(0),
    distance_(distance),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL),
    AMIPtr_(NULL),
    AMIReverse_(false),
    surfPtr_(NULL),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    sampleRegion_
    (
        dict.lookupOrDefault
        (
            "sampleRegion",
            patch_.boundaryMesh().mesh().name()
        )
    ),
    mode_(sampleModeNames_.read(dict.lookup("sampleMode"))),
    samplePatch_(dict.lookup("samplePatch")),
    offsetMode_(UNIFORM),
    offset_(vector::zero),
    offsets_(0),
    distance_(0.0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL),
    AMIPtr_(NULL),
    AMIReverse_(dict.lookupOrDefault<bool>("flipNormals", false)),
    surfPtr_(NULL),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (dict.found("offsetMode"))
    {
        offsetMode_ = offsetModeNames_.read(dict.lookup("offsetMode"));

        switch(offsetMode_)
        {
            case UNIFORM:
            {
                offset_ = point(dict.lookup("offset"));
            }
            break;

            case NONUNIFORM:
            {
                //offsets_ = pointField(dict.lookup("offsets"));
                offsets_ = readListOrField("offsets", dict, patch_.size());
            }
            break;

            case NORMAL:
            {
                distance_ = readScalar(dict.lookup("distance"));
            }
            break;
        }
    }
    else if (dict.found("offset"))
    {
        offsetMode_ = UNIFORM;
        offset_ = point(dict.lookup("offset"));
    }
    else if (dict.found("offsets"))
    {
        offsetMode_ = NONUNIFORM;
        //offsets_ = pointField(dict.lookup("offsets"));
        offsets_ = readListOrField("offsets", dict, patch_.size());
    }
    else
    {
        FatalIOErrorIn
        (
            "mappedPatchBase::mappedPatchBase\n"
            "(\n"
            "    const polyPatch&,\n"
            "    const dictionary&\n"
            ")\n",
            dict
        )   << "Please supply the offsetMode as one of "
            << NamedEnum<offsetMode, 3>::words()
            << exit(FatalIOError);
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb
)
:
    patch_(pp),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_(mpb.offsets_),
    distance_(mpb.distance_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(NULL),
    AMIPtr_(NULL),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(NULL),
    surfDict_(mpb.surfDict_)
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb,
    const labelUList& mapAddressing
)
:
    patch_(pp),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_(mpb.offsets_, mapAddressing),
    distance_(mpb.distance_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(NULL),
    AMIPtr_(NULL),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(NULL),
    surfDict_(mpb.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::~mappedPatchBase()
{
    clearOut();
}


void Foam::mappedPatchBase::clearOut()
{
    mapPtr_.clear();
    AMIPtr_.clear();
    surfPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedPatchBase::sampleMesh() const
{
    return patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
    (
        sampleRegion_
    );
}


const Foam::polyPatch& Foam::mappedPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    const label patchI = nbrMesh.boundaryMesh().findPatchID(samplePatch_);

    if (patchI == -1)
    {
        FatalErrorIn("mappedPatchBase::samplePolyPatch()")
            << "Cannot find patch " << samplePatch_
            << " in region " << sampleRegion_ << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchI];
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints() const
{
    tmp<pointField> tfld(new pointField(patch_.faceCentres()));
    pointField& fld = tfld();

    switch (offsetMode_)
    {
        case UNIFORM:
        {
            fld += offset_;
            break;
        }
        case NONUNIFORM:
        {
            fld += offsets_;
            break;
        }
        case NORMAL:
        {
            // Get outwards pointing normal
            vectorField n(patch_.faceAreas());
            n /= mag(n);

            fld += distance_*n;
            break;
        }
    }

    return tfld;
}


void Foam::mappedPatchBase::write(Ostream& os) const
{
    os.writeKeyword("sampleMode") << sampleModeNames_[mode_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("sampleRegion") << sampleRegion_
        << token::END_STATEMENT << nl;
    os.writeKeyword("samplePatch") << samplePatch_
        << token::END_STATEMENT << nl;

    os.writeKeyword("offsetMode") << offsetModeNames_[offsetMode_]
        << token::END_STATEMENT << nl;

    switch (offsetMode_)
    {
        case UNIFORM:
        {
            os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
            break;
        }
        case NONUNIFORM:
        {
            offsets_.writeEntry("offsets", os);
            break;
        }
        case NORMAL:
        {
            os.writeKeyword("distance") << distance_ << token::END_STATEMENT
                << nl;
            break;
        }
    }

    if (mode_ == NEARESTPATCHFACEAMI)
    {
        if (AMIReverse_)
        {
            os.writeKeyword("flipNormals") << AMIReverse_
                << token::END_STATEMENT << nl;
        }

        if (!surfDict_.empty())
        {
            os.writeKeyword(surfDict_.dictName());
            os  << surfDict_;
        }
    }
}


// ************************************************************************* //
