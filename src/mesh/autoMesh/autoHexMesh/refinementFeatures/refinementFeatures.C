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

#include "refinementFeatures.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementFeatures::refinementFeatures
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts
)
:
    PtrList<featureEdgeMesh>(featDicts.size()),
    levels_(featDicts.size()),
    edgeTrees_(featDicts.size()),
    pointTrees_(featDicts.size())
{
    // Read features

    forAll(featDicts, i)
    {
        const dictionary& dict = featDicts[i];

        fileName featFileName(dict.lookup("file"));

        set
        (
            i,
            new featureEdgeMesh
            (
                IOobject
                (
                    featFileName,                       // name
                    io.time().constant(),               // instance
                    "triSurface",                       // local
                    io.time(),                          // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );

        const featureEdgeMesh& eMesh = operator[](i);

        //eMesh.mergePoints(meshRefiner_.mergeDistance());
        levels_[i] = readLabel(dict.lookup("level"));

        Info<< "Refinement level " << levels_[i]
            << " for all cells crossed by feature " << featFileName
            << " (" << eMesh.points().size() << " points, "
            << eMesh.edges().size() << " edges)." << endl;
    }


    // Search engines

    forAll(*this, i)
    {
        const featureEdgeMesh& eMesh = operator[](i);
        const pointField& points = eMesh.points();
        const edgeList& edges = eMesh.edges();

        // Calculate bb of all points
        treeBoundBox bb(points);

        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        bb = bb.extend(rndGen, 1E-4);
        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        edgeTrees_.set
        (
            i,
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,                  // do not cache bb
                    edges,
                    points,
                    identity(edges.size())
                ),
                bb,     // overall search domain
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );


        // Detect feature points from edges.
        const labelListList& pointEdges = eMesh.pointEdges();
        DynamicList<label> featurePoints;
        forAll(pointEdges, pointI)
        {
            if (pointEdges[pointI].size() > 2)
            {
                featurePoints.append(pointI);
            }
        }

        Info<< "Detected " << featurePoints.size()
            << " featurePoints out of " << points.size()
            << " on feature " << eMesh.name() << endl;

        pointTrees_.set
        (
            i,
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint(points, featurePoints),
                bb,     // overall search domain
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::refinementFeatures::findNearestEdge
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());

    forAll(edgeTrees_, featI)
    {
        const indexedOctree<treeDataEdge>& tree = edgeTrees_[featI];

        if (tree.shapes().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearInfo[sampleI].hit())
                {
                    distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearInfo[sampleI] = info;
                    nearFeature[sampleI] = featI;
                }
            }
        }
    }
}


void Foam::refinementFeatures::findNearestPoint
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    labelList& nearIndex
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearIndex.setSize(samples.size());
    nearIndex = -1;

    forAll(pointTrees_, featI)
    {
        const indexedOctree<treeDataPoint>& tree = pointTrees_[featI];

        if (tree.shapes().pointLabels().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearFeature[sampleI] != -1)
                {
                    label nearFeatI = nearFeature[sampleI];
                    const indexedOctree<treeDataPoint>& nearTree =
                        pointTrees_[nearFeatI];
                    label featPointI =
                        nearTree.shapes().pointLabels()[nearIndex[sampleI]];
                    const point& featPt =
                        operator[](nearFeatI).points()[featPointI];
                    distSqr = magSqr(featPt-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[sampleI] = featI;
                    nearIndex[sampleI] = info.index();
                }
            }
        }
    }
}


// ************************************************************************* //
