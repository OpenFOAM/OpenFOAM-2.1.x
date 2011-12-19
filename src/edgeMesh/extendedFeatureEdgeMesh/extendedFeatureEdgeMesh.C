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

#include "extendedFeatureEdgeMesh.H"
#include "triSurface.H"
#include "Random.H"
#include "Time.H"
#include "meshTools.H"
#include "linePointRef.H"
#include "ListListOps.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::extendedFeatureEdgeMesh, 0);

Foam::scalar Foam::extendedFeatureEdgeMesh::cosNormalAngleTol_ =
    Foam::cos(degToRad(0.1));


Foam::label Foam::extendedFeatureEdgeMesh::convexStart_ = 0;


Foam::label Foam::extendedFeatureEdgeMesh::externalStart_ = 0;


Foam::label Foam::extendedFeatureEdgeMesh::nPointTypes = 4;


Foam::label Foam::extendedFeatureEdgeMesh::nEdgeTypes = 5;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh(const IOobject& io)
:
    regIOobject(io),
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    edgeDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    regionEdges_(0),
    edgeTree_(),
    edgeTreesByType_()
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningIn
            (
                "extendedFeatureEdgeMesh::extendedFeatureEdgeMesh"
                "(const IOobject&)"
            )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }

        Istream& is = readStream(typeName);

        is  >> *this
            >> concaveStart_
            >> mixedStart_
            >> nonFeatureStart_
            >> internalStart_
            >> flatStart_
            >> openStart_
            >> multipleStart_
            >> normals_
            >> edgeNormals_
            >> featurePointNormals_
            >> regionEdges_;

        close();

        {
            // Calculate edgeDirections

            const edgeList& eds(edges());

            const pointField& pts(points());

            edgeDirections_.setSize(eds.size());

            forAll(eds, eI)
            {
                edgeDirections_[eI] = eds[eI].vec(pts);
            }

            edgeDirections_ /= mag(edgeDirections_);
        }
    }

    if (debug)
    {
        Pout<< "extendedFeatureEdgeMesh::extendedFeatureEdgeMesh :"
            << " constructed from IOobject :"
            << " points:" << points().size()
            << " edges:" << edges().size()
            << endl;
    }
}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const extendedFeatureEdgeMesh& fem
)
:
    regIOobject(io),
    edgeMesh(fem),
    concaveStart_(fem.concaveStart()),
    mixedStart_(fem.mixedStart()),
    nonFeatureStart_(fem.nonFeatureStart()),
    internalStart_(fem.internalStart()),
    flatStart_(fem.flatStart()),
    openStart_(fem.openStart()),
    multipleStart_(fem.multipleStart()),
    normals_(fem.normals()),
    edgeDirections_(fem.edgeDirections()),
    edgeNormals_(fem.edgeNormals()),
    featurePointNormals_(fem.featurePointNormals()),
    regionEdges_(fem.regionEdges()),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const Xfer<pointField>& pointLst,
    const Xfer<edgeList>& edgeLst
)
:
    regIOobject(io),
    edgeMesh(pointLst, edgeLst),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    edgeDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    regionEdges_(0),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const surfaceFeatures& sFeat,
    const objectRegistry& obr,
    const fileName& sFeatFileName
)
:
    regIOobject
    (
        IOobject
        (
            sFeatFileName,
            obr.time().constant(),
            "extendedFeatureEdgeMesh",
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(-1),
    mixedStart_(-1),
    nonFeatureStart_(-1),
    internalStart_(-1),
    flatStart_(-1),
    openStart_(-1),
    multipleStart_(-1),
    normals_(0),
    edgeDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    regionEdges_(0),
    edgeTree_(),
    edgeTreesByType_()
{
    // Extract and reorder the data from surfaceFeatures

    // References to the surfaceFeatures data
    const triSurface& surf(sFeat.surface());
    const pointField& sFeatLocalPts(surf.localPoints());
    const edgeList& sFeatEds(surf.edges());

    // Filling the extendedFeatureEdgeMesh with the raw geometrical data.

    label nFeatEds = sFeat.featureEdges().size();

    DynamicList<point> tmpPts;
    edgeList eds(nFeatEds);
    DynamicList<vector> norms;
    vectorField edgeDirections(nFeatEds);
    labelListList edgeNormals(nFeatEds);
    DynamicList<label> regionEdges;

    // Mapping between old and new indices, there is entry in the map for each
    // of sFeatLocalPts, -1 means that this point hasn't been used (yet), >= 0
    // corresponds to the index
    labelList pointMap(sFeatLocalPts.size(), -1);

    // Noting when the normal of a face has been used so not to duplicate
    labelList faceMap(surf.size(), -1);

    // Collecting the status of edge for subsequent sorting
    List<edgeStatus> edStatus(nFeatEds, NONE);

    forAll(sFeat.featurePoints(), i)
    {
        label sFPI = sFeat.featurePoints()[i];

        tmpPts.append(sFeatLocalPts[sFPI]);

        pointMap[sFPI] = tmpPts.size() - 1;
    }

    // All feature points have been added
    nonFeatureStart_ = tmpPts.size();

    forAll(sFeat.featureEdges(), i)
    {
        label sFEI = sFeat.featureEdges()[i];

        const edge& fE(sFeatEds[sFEI]);

        // Check to see if the points have been already used
        if (pointMap[fE.start()] == -1)
        {
             tmpPts.append(sFeatLocalPts[fE.start()]);

             pointMap[fE.start()] = tmpPts.size() - 1;
        }

        eds[i].start() = pointMap[fE.start()];

        if (pointMap[fE.end()] == -1)
        {
            tmpPts.append(sFeatLocalPts[fE.end()]);

            pointMap[fE.end()] = tmpPts.size() - 1;
        }

        eds[i].end() = pointMap[fE.end()];

        // Pick up the faces adjacent to the feature edge
        const labelList& eFaces = surf.edgeFaces()[sFEI];

        edgeNormals[i].setSize(eFaces.size());

        forAll(eFaces, j)
        {
            label eFI = eFaces[j];

            // Check to see if the points have been already used
            if (faceMap[eFI] == -1)
            {
                norms.append(surf.faceNormals()[eFI]);

                faceMap[eFI] = norms.size() - 1;
            }

            edgeNormals[i][j] = faceMap[eFI];
        }

        vector fC0tofC1(vector::zero);

        if (eFaces.size() == 2)
        {
            fC0tofC1 =
                surf[eFaces[1]].centre(surf.points())
              - surf[eFaces[0]].centre(surf.points());
        }

        edStatus[i] = classifyEdge(norms, edgeNormals[i], fC0tofC1);

        edgeDirections[i] = fE.vec(sFeatLocalPts);

        if (i < sFeat.nRegionEdges())
        {
            regionEdges.append(i);
        }
    }

    // Reorder the edges by classification

    List<DynamicList<label> > allEds(nEdgeTypes);

    DynamicList<label>& externalEds(allEds[0]);
    DynamicList<label>& internalEds(allEds[1]);
    DynamicList<label>& flatEds(allEds[2]);
    DynamicList<label>& openEds(allEds[3]);
    DynamicList<label>& multipleEds(allEds[4]);

    forAll(eds, i)
    {
        edgeStatus eStat = edStatus[i];

        if (eStat == EXTERNAL)
        {
            externalEds.append(i);
        }
        else if (eStat == INTERNAL)
        {
            internalEds.append(i);
        }
        else if (eStat == FLAT)
        {
            flatEds.append(i);
        }
        else if (eStat == OPEN)
        {
            openEds.append(i);
        }
        else if (eStat == MULTIPLE)
        {
            multipleEds.append(i);
        }
        else if (eStat == NONE)
        {
            FatalErrorIn
            (
                "Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh"
                "("
                "    const surfaceFeatures& sFeat,"
                "    const objectRegistry& obr,"
                "    const fileName& sFeatFileName"
                ")"
            )
                << nl << "classifyEdge returned NONE on edge "
                << eds[i]
                << ". There is a problem with definition of this edge."
                << nl << abort(FatalError);
        }
    }

    internalStart_ = externalEds.size();
    flatStart_ = internalStart_ + internalEds.size();
    openStart_ = flatStart_ + flatEds.size();
    multipleStart_ = openStart_ + openEds.size();

    labelList edMap
    (
        ListListOps::combine<labelList>
        (
            allEds,
            accessOp<labelList>()
        )
    );

    edMap = invert(edMap.size(), edMap);

    inplaceReorder(edMap, eds);
    inplaceReorder(edMap, edStatus);
    inplaceReorder(edMap, edgeDirections);
    inplaceReorder(edMap, edgeNormals);
    inplaceRenumber(edMap, regionEdges);

    pointField pts(tmpPts);

    // Initialise the edgeMesh
    edgeMesh::operator=(edgeMesh(pts, eds));

    // Initialise sorted edge related data
    edgeDirections_ = edgeDirections/mag(edgeDirections);
    edgeNormals_ = edgeNormals;
    regionEdges_ = regionEdges;

    // Normals are all now found and indirectly addressed, can also be stored
    normals_ = vectorField(norms);

    // Reorder the feature points by classification

    List<DynamicList<label> > allPts(3);

    DynamicList<label>& convexPts(allPts[0]);
    DynamicList<label>& concavePts(allPts[1]);
    DynamicList<label>& mixedPts(allPts[2]);

    for (label i = 0; i < nonFeatureStart_; i++)
    {
        pointStatus ptStatus = classifyFeaturePoint(i);

        if (ptStatus == CONVEX)
        {
            convexPts.append(i);
        }
        else if (ptStatus == CONCAVE)
        {
            concavePts.append(i);
        }
        else if (ptStatus == MIXED)
        {
            mixedPts.append(i);
        }
        else if (ptStatus == NONFEATURE)
        {
            FatalErrorIn
            (
                "Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh"
                "("
                "    const surfaceFeatures& sFeat,"
                "    const objectRegistry& obr,"
                "    const fileName& sFeatFileName"
                ")"
            )
                << nl << "classifyFeaturePoint returned NONFEATURE on point at "
                << points()[i]
                << ". There is a problem with definition of this feature point."
                << nl << abort(FatalError);
        }
    }

    concaveStart_ = convexPts.size();
    mixedStart_ = concaveStart_ + concavePts.size();

    labelList ftPtMap
    (
        ListListOps::combine<labelList>
        (
            allPts,
            accessOp<labelList>()
        )
    );

    ftPtMap = invert(ftPtMap.size(), ftPtMap);

    // Creating the ptMap from the ftPtMap with identity values up to the size
    // of pts to create an oldToNew map for inplaceReorder

    labelList ptMap(identity(pts.size()));

    forAll(ftPtMap, i)
    {
        ptMap[i] = ftPtMap[i];
    }

    inplaceReorder(ptMap, pts);

    forAll(eds, i)
    {
        inplaceRenumber(ptMap, eds[i]);
    }

    // Reinitialise the edgeMesh with sorted feature points and
    // renumbered edges
    edgeMesh::operator=(edgeMesh(pts, eds));

    // Generate the featurePointNormals

    labelListList featurePointNormals(nonFeatureStart_);

    for (label i = 0; i < nonFeatureStart_; i++)
    {
        DynamicList<label> tmpFtPtNorms;

        const labelList& ptEds = pointEdges()[i];

        forAll(ptEds, j)
        {
            const labelList& ptEdNorms(edgeNormals[ptEds[j]]);

            forAll(ptEdNorms, k)
            {
                if (findIndex(tmpFtPtNorms, ptEdNorms[k]) == -1)
                {
                    bool addNormal = true;

                    // Check that the normal direction is unique at this feature
                    forAll(tmpFtPtNorms, q)
                    {
                        if
                        (
                            (normals_[ptEdNorms[k]] & normals_[tmpFtPtNorms[q]])
                          > cosNormalAngleTol_
                        )
                        {
                            // Parallel to an existing normal, do not add
                            addNormal = false;

                            break;
                        }
                    }

                    if (addNormal)
                    {
                        tmpFtPtNorms.append(ptEdNorms[k]);
                    }
                }
            }
        }

        featurePointNormals[i] = tmpFtPtNorms;
    }

    featurePointNormals_ = featurePointNormals;
}


Foam::extendedFeatureEdgeMesh::extendedFeatureEdgeMesh
(
    const IOobject& io,
    const pointField& pts,
    const edgeList& eds,
    label concaveStart,
    label mixedStart,
    label nonFeatureStart,
    label internalStart,
    label flatStart,
    label openStart,
    label multipleStart,
    const vectorField& normals,
    const vectorField& edgeDirections,
    const labelListList& edgeNormals,
    const labelListList& featurePointNormals,
    const labelList& regionEdges
)
:
    regIOobject(io),
    edgeMesh(pts, eds),
    concaveStart_(concaveStart),
    mixedStart_(mixedStart),
    nonFeatureStart_(nonFeatureStart),
    internalStart_(internalStart),
    flatStart_(flatStart),
    openStart_(openStart),
    multipleStart_(multipleStart),
    normals_(normals),
    edgeDirections_(edgeDirections),
    edgeNormals_(edgeNormals),
    featurePointNormals_(featurePointNormals),
    regionEdges_(regionEdges),
    edgeTree_(),
    edgeTreesByType_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::~extendedFeatureEdgeMesh()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::extendedFeatureEdgeMesh::pointStatus
Foam::extendedFeatureEdgeMesh::classifyFeaturePoint
(
    label ptI
) const
{
    labelList ptEds(pointEdges()[ptI]);

    label nPtEds = ptEds.size();
    label nExternal = 0;
    label nInternal = 0;

    if (nPtEds == 0)
    {
        // There are no edges attached to the point, this is a problem
        return NONFEATURE;
    }

    forAll(ptEds, i)
    {
        edgeStatus edStat = getEdgeStatus(ptEds[i]);

        if (edStat == EXTERNAL)
        {
            nExternal++;
        }
        else if (edStat == INTERNAL)
        {
            nInternal++;
        }
    }

    if (nExternal == nPtEds)
    {
        return CONVEX;
    }
    else if (nInternal == nPtEds)
    {
        return CONCAVE;
    }
    else
    {
        return MIXED;
    }
}


Foam::extendedFeatureEdgeMesh::edgeStatus
Foam::extendedFeatureEdgeMesh::classifyEdge
(
    const List<vector>& norms,
    const labelList& edNorms,
    const vector& fC0tofC1
) const
{
    label nEdNorms = edNorms.size();

    if (nEdNorms == 1)
    {
        return OPEN;
    }
    else if (nEdNorms == 2)
    {
        const vector n0(norms[edNorms[0]]);
        const vector n1(norms[edNorms[1]]);

        if ((n0 & n1) > cosNormalAngleTol_)
        {
            return FLAT;
        }
        else if ((fC0tofC1 & n0) > 0.0)
        {
            return INTERNAL;
        }
        else
        {
            return EXTERNAL;
        }
    }
    else if (nEdNorms > 2)
    {
        return MULTIPLE;
    }
    else
    {
        // There is a problem - the edge has no normals
        return NONE;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::extendedFeatureEdgeMesh::nearestFeatureEdge
(
    const point& sample,
    scalar searchDistSqr,
    pointIndexHit& info
) const
{
    info = edgeTree().findNearest
    (
        sample,
        searchDistSqr
    );
}


void Foam::extendedFeatureEdgeMesh::nearestFeatureEdge
(
    const pointField& samples,
    const scalarField& searchDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        nearestFeatureEdge
        (
            samples[i],
            searchDistSqr[i],
            info[i]
        );
    }
}


void Foam::extendedFeatureEdgeMesh::nearestFeatureEdgeByType
(
    const point& sample,
    const scalarField& searchDistSqr,
    List<pointIndexHit>& info
) const
{
    const PtrList<indexedOctree<treeDataEdge> >& edgeTrees = edgeTreesByType();

    info.setSize(edgeTrees.size());

    labelList sliceStarts(edgeTrees.size());

    sliceStarts[0] = externalStart_;
    sliceStarts[1] = internalStart_;
    sliceStarts[2] = flatStart_;
    sliceStarts[3] = openStart_;
    sliceStarts[4] = multipleStart_;

    forAll(edgeTrees, i)
    {
        info[i] = edgeTrees[i].findNearest
        (
            sample,
            searchDistSqr[i]
        );

        // The index returned by the indexedOctree is local to the slice of
        // edges it was supplied with, return the index to the value in the
        // complete edge list

        info[i].setIndex(info[i].index() + sliceStarts[i]);
    }
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::extendedFeatureEdgeMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1E-4)
        );

        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        labelList allEdges(identity(edges().size()));

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    points(),       // points
                    allEdges        // selected edges
                ),
                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return edgeTree_();
}


const Foam::PtrList<Foam::indexedOctree<Foam::treeDataEdge> >&
Foam::extendedFeatureEdgeMesh::edgeTreesByType() const
{
    if (edgeTreesByType_.size() == 0)
    {
        edgeTreesByType_.setSize(nEdgeTypes);

        Random rndGen(872141);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1E-4)
        );

        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        labelListList sliceEdges(nEdgeTypes);

        // External edges
        sliceEdges[0] =
            identity(internalStart_ - externalStart_) + externalStart_;

        // Internal edges
        sliceEdges[1] = identity(flatStart_ - internalStart_) + internalStart_;

        // Flat edges
        sliceEdges[2] = identity(openStart_ - flatStart_) + flatStart_;

        // Open edges
        sliceEdges[3] = identity(multipleStart_ - openStart_) + openStart_;

        // Multiple edges
        sliceEdges[4] =
            identity(edges().size() - multipleStart_) + multipleStart_;

        forAll(edgeTreesByType_, i)
        {
            edgeTreesByType_.set
            (
                i,
                new indexedOctree<treeDataEdge>
                (
                    treeDataEdge
                    (
                        false,          // cachebb
                        edges(),        // edges
                        points(),       // points
                        sliceEdges[i]   // selected edges
                    ),
                    bb,     // bb
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                )
            );
        }
    }

    return edgeTreesByType_;
}


void Foam::extendedFeatureEdgeMesh::writeObj
(
    const fileName& prefix
) const
{
    Pout<< nl << "Writing extendedFeatureEdgeMesh components to " << prefix
        << endl;

    label verti = 0;

    edgeMesh::write(prefix + "_edgeMesh.obj");

    OFstream convexFtPtStr(prefix + "_convexFeaturePts.obj");
    Pout<< "Writing convex feature points to " << convexFtPtStr.name() << endl;

    for(label i = 0; i < concaveStart_; i++)
    {
        meshTools::writeOBJ(convexFtPtStr, points()[i]);
    }

    OFstream concaveFtPtStr(prefix + "_concaveFeaturePts.obj");
    Pout<< "Writing concave feature points to "
        << concaveFtPtStr.name() << endl;

    for(label i = concaveStart_; i < mixedStart_; i++)
    {
        meshTools::writeOBJ(concaveFtPtStr, points()[i]);
    }

    OFstream mixedFtPtStr(prefix + "_mixedFeaturePts.obj");
    Pout<< "Writing mixed feature points to " << mixedFtPtStr.name() << endl;

    for(label i = mixedStart_; i < nonFeatureStart_; i++)
    {
        meshTools::writeOBJ(mixedFtPtStr, points()[i]);
    }

    OFstream mixedFtPtStructureStr(prefix + "_mixedFeaturePtsStructure.obj");
    Pout<< "Writing mixed feature point structure to "
        << mixedFtPtStructureStr.name() << endl;

    verti = 0;
    for(label i = mixedStart_; i < nonFeatureStart_; i++)
    {
        const labelList& ptEds = pointEdges()[i];

        forAll(ptEds, j)
        {
            const edge& e = edges()[ptEds[j]];

            meshTools::writeOBJ(mixedFtPtStructureStr, points()[e[0]]); verti++;
            meshTools::writeOBJ(mixedFtPtStructureStr, points()[e[1]]); verti++;
            mixedFtPtStructureStr << "l " << verti-1 << ' ' << verti << endl;
        }
    }

    OFstream externalStr(prefix + "_externalEdges.obj");
    Pout<< "Writing external edges to " << externalStr.name() << endl;

    verti = 0;
    for (label i = externalStart_; i < internalStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(externalStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(externalStr, points()[e[1]]); verti++;
        externalStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream internalStr(prefix + "_internalEdges.obj");
    Pout<< "Writing internal edges to " << internalStr.name() << endl;

    verti = 0;
    for (label i = internalStart_; i < flatStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(internalStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(internalStr, points()[e[1]]); verti++;
        internalStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream flatStr(prefix + "_flatEdges.obj");
    Pout<< "Writing flat edges to " << flatStr.name() << endl;

    verti = 0;
    for (label i = flatStart_; i < openStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(flatStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(flatStr, points()[e[1]]); verti++;
        flatStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream openStr(prefix + "_openEdges.obj");
    Pout<< "Writing open edges to " << openStr.name() << endl;

    verti = 0;
    for (label i = openStart_; i < multipleStart_; i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(openStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(openStr, points()[e[1]]); verti++;
        openStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream multipleStr(prefix + "_multipleEdges.obj");
    Pout<< "Writing multiple edges to " << multipleStr.name() << endl;

    verti = 0;
    for (label i = multipleStart_; i < edges().size(); i++)
    {
        const edge& e = edges()[i];

        meshTools::writeOBJ(multipleStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(multipleStr, points()[e[1]]); verti++;
        multipleStr << "l " << verti-1 << ' ' << verti << endl;
    }

    OFstream regionStr(prefix + "_regionEdges.obj");
    Pout<< "Writing region edges to " << regionStr.name() << endl;

    verti = 0;
    forAll(regionEdges_, i)
    {
        const edge& e = edges()[regionEdges_[i]];

        meshTools::writeOBJ(regionStr, points()[e[0]]); verti++;
        meshTools::writeOBJ(regionStr, points()[e[1]]); verti++;
        regionStr << "l " << verti-1 << ' ' << verti << endl;
    }
}


bool Foam::extendedFeatureEdgeMesh::writeData(Ostream& os) const
{
    os  << "// points, edges, concaveStart, mixedStart, nonFeatureStart, " << nl
        << "// internalStart, flatStart, openStart, multipleStart, " << nl
        << "// normals, edgeNormals, featurePointNormals, regionEdges" << nl
        << endl;

    os  << points() << nl
        << edges() << nl
        << concaveStart_ << token::SPACE
        << mixedStart_ << token::SPACE
        << nonFeatureStart_ << token::SPACE
        << internalStart_ << token::SPACE
        << flatStart_ << token::SPACE
        << openStart_ << token::SPACE
        << multipleStart_ << nl
        << normals_ << nl
        << edgeNormals_ << nl
        << featurePointNormals_ << nl
        << regionEdges_
        << endl;

    return os.good();
}


// ************************************************************************* //
