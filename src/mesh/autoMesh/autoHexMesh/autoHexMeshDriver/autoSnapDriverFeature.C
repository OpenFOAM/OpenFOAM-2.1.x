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

\*----------------------------------------------------------------------------*/

#include "autoSnapDriver.H"
#include "polyTopoChange.H"
#include "OFstream.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "motionSmoother.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "unitConversion.H"
#include "plane.H"
#include "featureEdgeMesh.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    class listTransform
    {
    public:

        void operator()
        (
            const vectorTensorTransform& vt,
            const bool forward,
            List<List<point> >& fld
        ) const
        {
            const tensor T
            (
                forward
              ? vt.R()
              : vt.R().T()
            );

            forAll(fld, i)
            {
                List<point>& elems = fld[i];
                forAll(elems, elemI)
                {
                    elems[elemI] = transform(T, elems[elemI]);
                }
            }
        }
    };

    class listPlusEqOp
    {
    public:

        void operator()
        (
            List<point>& x,
            const List<point>& y
        ) const
        {
            label sz = x.size();
            x.setSize(sz+y.size());
            forAll(y, i)
            {
                x[sz++] = y[i];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::autoSnapDriver::smoothAndConstrain
(
    const indirectPrimitivePatch& pp,
    const List<pointConstraint>& constraints,
    vectorField& disp
) const
{
    for (label avgIter = 0; avgIter < 20; avgIter++)
    {
        // Calculate average displacement of neighbours
        vectorField dispAvg(pp.nPoints(), vector::zero);

        const labelListList& pointEdges = pp.pointEdges();

        forAll(pointEdges, pointI)
        {
            const labelList& pEdges = pointEdges[pointI];

            forAll(pEdges, i)
            {
                const edge& e = pp.edges()[pEdges[i]];
                label nbrPointI = e.otherVertex(pointI);
                dispAvg[pointI] += disp[nbrPointI];
            }
            dispAvg[pointI] /= pEdges.size();
        }

        // Constraints
        forAll(constraints, pointI)
        {
            if (constraints[pointI].first() == 0)
            {
                // Mix my displacement with neighbours' displacement
                disp[pointI] = 0.5*disp[pointI] + 0.5*dispAvg[pointI];
            }
        }
    }
}


void Foam::autoSnapDriver::calcNearestFace
(
    const label iter,
    const indirectPrimitivePatch& pp,
    vectorField& faceDisp,
    vectorField& faceSurfaceNormal,
    vectorField& faceRotation
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Displacement and orientation per pp face.
    faceDisp.setSize(pp.size());
    faceDisp = vector::zero;
    faceSurfaceNormal.setSize(pp.size());
    faceSurfaceNormal = vector::zero;

    // Divide surfaces into zoned and unzoned
    labelList zonedSurfaces = surfaces.getNamedSurfaces();
    labelList unzonedSurfaces = surfaces.getUnnamedSurfaces();

    // Per pp face the current surface snapped to
    labelList snapSurf(pp.size(), -1);


    // Do zoned surfaces
    // ~~~~~~~~~~~~~~~~~
    // Zoned faces only attract to corresponding surface

    // Extract faces per zone
    const wordList& faceZoneNames = surfaces.faceZoneNames();

    forAll(zonedSurfaces, i)
    {
        label zoneSurfI = zonedSurfaces[i];

        // Get indices of faces on pp that are also in zone
        label zoneI = mesh.faceZones().findZoneID(faceZoneNames[zoneSurfI]);
        if (zoneI == -1)
        {
            FatalErrorIn
            (
                "autoSnapDriver::calcNearestFace(..)"
            )   << "Problem. Cannot find zone " << faceZoneNames[zoneSurfI]
                << exit(FatalError);
        }
        const faceZone& fZone = mesh.faceZones()[zoneI];
        PackedBoolList isZonedFace(mesh.nFaces());
        forAll(fZone, i)
        {
            isZonedFace[fZone[i]] = 1;
        }

        DynamicList<label> ppFaces(fZone.size());
        DynamicList<label> meshFaces(fZone.size());
        forAll(pp.addressing(), i)
        {
            if (isZonedFace[pp.addressing()[i]])
            {
                snapSurf[i] = zoneSurfI;
                ppFaces.append(i);
                meshFaces.append(pp.addressing()[i]);
            }
        }

        //Pout<< "For faceZone " << fZone.name()
        //    << " found " << ppFaces.size() << " out of " << pp.size()
        //    << endl;

        pointField fc
        (
            indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), meshFaces),
                mesh.points()
            ).faceCentres()
        );

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        labelList hitRegion;
        vectorField hitNormal;
        surfaces.findNearestRegion
        (
            labelList(1, zoneSurfI),
            fc,
            sqr(scalarField(fc.size(), GREAT)),// sqr of attract dist
            hitSurface,
            hitInfo,
            hitRegion,
            hitNormal
        );

        forAll(hitInfo, hitI)
        {
            if (hitInfo[hitI].hit())
            {
                label faceI = ppFaces[hitI];
                faceDisp[faceI] = hitInfo[hitI].hitPoint() - fc[hitI];
                faceSurfaceNormal[faceI] = hitNormal[hitI];
            }
        }
    }


    // Do unzoned surfaces
    // ~~~~~~~~~~~~~~~~~~~
    // Unzoned faces attract to any unzoned surface

    DynamicList<label> ppFaces(pp.size());
    DynamicList<label> meshFaces(pp.size());
    forAll(pp.addressing(), i)
    {
        if (snapSurf[i] == -1)
        {
            ppFaces.append(i);
            meshFaces.append(pp.addressing()[i]);
        }
    }
    //Pout<< "Found " << ppFaces.size() << " unzoned faces out of " << pp.size()
    //    << endl;

    pointField fc
    (
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), meshFaces),
            mesh.points()
        ).faceCentres()
    );

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    labelList hitRegion;
    vectorField hitNormal;
    surfaces.findNearestRegion
    (
        unzonedSurfaces,
        fc,
        sqr(scalarField(fc.size(), GREAT)),// sqr of attract dist
        hitSurface,
        hitInfo,
        hitRegion,
        hitNormal
    );

    forAll(hitInfo, hitI)
    {
        if (hitInfo[hitI].hit())
        {
            label faceI = ppFaces[hitI];
            faceDisp[faceI] = hitInfo[hitI].hitPoint() - fc[hitI];
            faceSurfaceNormal[faceI] = hitNormal[hitI];
        }
    }


    // Determine rotation
    // ~~~~~~~~~~~~~~~~~~

    // Determine rotation axis
    faceRotation.setSize(pp.size());
    faceRotation = vector::zero;

    forAll(faceRotation, faceI)
    {
        // Note: extend to >180 degrees checking
        faceRotation[faceI] =
            pp.faceNormals()[faceI]
          ^ -faceSurfaceNormal[faceI];
    }

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        dumpMove
        (
            mesh.time().path()
          / "faceDisp_" + name(iter) + ".obj",
            pp.faceCentres(),
            pp.faceCentres() + faceDisp
        );
        dumpMove
        (
            mesh.time().path()
          / "faceRotation_" + name(iter) + ".obj",
            pp.faceCentres(),
            pp.faceCentres() + faceRotation
        );
    }
}


// Gets passed in offset to nearest point on feature edge. Calculates
// if the point has a different number of faces on either side of the feature
// and if so attracts the point to that non-dominant plane.
void Foam::autoSnapDriver::correctAttraction
(
    const DynamicList<point>& surfacePoints,
    const DynamicList<label>& surfaceCount,
    const point& edgePt,
    const vector& edgeNormal,       // normalised normal
    const point& pt,

    vector& edgeOffset              // offset from pt to point on edge
) const
{
    // Tangential component along edge
    scalar tang = ((pt-edgePt)&edgeNormal);

    labelList order;
    Foam::sortedOrder(surfaceCount, order);

    if (order[0] < order[1])
    {
        // There is a non-dominant plane. Use the point on the plane to
        // attract to.
        vector attractD = surfacePoints[order[0]]-edgePt;
        // Tangential component along edge
        scalar tang2 = (attractD&edgeNormal);
        // Normal component
        attractD -= tang2*edgeNormal;
        // Calculate fraction of normal distances
        scalar magAttractD = mag(attractD);
        scalar fraction = magAttractD/(magAttractD+mag(edgeOffset));

        point linePt =
            edgePt
          + ((1.0-fraction)*tang2 + fraction*tang)*edgeNormal;
        edgeOffset = linePt-pt;
    }
}


void Foam::autoSnapDriver::binFeatureFace
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalar snapDist,

    const point& fc,
    const vector& faceSurfaceNormal,
    const vector& faceDisp,

    DynamicList<point>& surfacePoints,
    DynamicList<vector>& surfaceNormals,
    DynamicList<label>& surfaceCount
) const
{
    // What to do with very far attraction? For now just ignore the face
    if (magSqr(faceDisp) < sqr(snapDist))
    {
        const point pt = fc + faceDisp;

        bool same = false;
        forAll(surfaceNormals, j)
        {
            scalar cosAngle = (faceSurfaceNormal&surfaceNormals[j]);

            if
            (
                (cosAngle >= featureCos)
             || (cosAngle < (-1+0.001)) // triangle baffles
            )
            {
                same = true;
                surfaceCount[j]++;
                break;
            }
        }

        if (!same)
        {
            // Now check if the planes go through the same edge or point

            if (surfacePoints.size() <= 1)
            {
                surfacePoints.append(pt);
                surfaceNormals.append(faceSurfaceNormal);
                surfaceCount.append(1);
            }
            else if (surfacePoints.size() == 2)
            {
                plane pl0(surfacePoints[0], surfaceNormals[0]);
                plane pl1(surfacePoints[1], surfaceNormals[1]);
                plane::ray r(pl0.planeIntersect(pl1));
                vector featureNormal = r.dir() / mag(r.dir());

                if (mag(faceSurfaceNormal&featureNormal) >= 0.001)
                {
                    // Definitely makes a feature point
                    surfacePoints.append(pt);
                    surfaceNormals.append(faceSurfaceNormal);
                    surfaceCount.append(1);
                }
            }
            else if (surfacePoints.size() == 3)
            {
                // Have already feature point. See if this new plane is the
                // same point or not.
                plane pl0(surfacePoints[0], surfaceNormals[0]);
                plane pl1(surfacePoints[1], surfaceNormals[1]);
                plane pl2(surfacePoints[2], surfaceNormals[2]);
                point p012(pl0.planePlaneIntersect(pl1, pl2));

                plane::ray r(pl0.planeIntersect(pl1));
                vector featureNormal = r.dir() / mag(r.dir());
                if (mag(faceSurfaceNormal&featureNormal) >= 0.001)
                {
                    plane pl3(pt, faceSurfaceNormal);
                    point p013(pl0.planePlaneIntersect(pl1, pl3));

                    if (mag(p012-p013) > 0.0001)    //TBD
                    {
                        // Different feature point
                        //Pout<< "** differing feature point :" << p012
                        //    << " and " << p013 << endl;
                        // Mark point as illegal (for now by setting size
                        // to 4)
                        surfacePoints.append(pt);
                        surfaceNormals.append(faceSurfaceNormal);
                        surfaceCount.append(1);
                    }
                }
            }
        }
    }
    //else
    //{
    //    Pout<< "binFeatureFace : for point:" << pp.localPoints()[pointI]
    //        << " found face-nearest point:"
    //        << fc + faceDisp
    //        << " at distance:" << mag(faceDisp)
    //        << " further away than:" << snapDist
    //        << endl;
    //}
}


// Check the faces surrounding a point. Bin them according to normal.
void Foam::autoSnapDriver::binFeatureFaces
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const label pointI,

    const List<List<point> >& pointFaceNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,

    DynamicList<point>& surfacePoints,
    DynamicList<vector>& surfaceNormals,
    DynamicList<label>& surfaceCount
) const
{
    const List<point>& pfNormals = pointFaceNormals[pointI];
    const List<point>& pfDisp = pointFaceDisp[pointI];
    const List<point>& pfCentres = pointFaceCentres[pointI];

    // Collect all different directions
    forAll(pfNormals, i)
    {
        binFeatureFace
        (
            iter,
            featureCos,

            pp,
            snapDist[pointI],

            pfCentres[i],
            pfNormals[i],
            pfDisp[i],

            surfacePoints,
            surfaceNormals,
            surfaceCount
        );
    }
}


// Special version that calculates attraction in one go
void Foam::autoSnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const List<List<point> >& pointFaceNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OFstream> featureEdgeStr;
    label featureEdgeVertI = 0;
    autoPtr<OFstream> featurePointStr;
    label featurePointVertI = 0;

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        featureEdgeStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featureEdge_" + name(iter) + ".obj"
            )
        );
        Pout<< nl << "Dumping feature-edge direction to "
            << featureEdgeStr().name() << endl;

        featurePointStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featurePoint_" + name(iter) + ".obj"
            )
        );
        Pout<< nl << "Dumping feature-point direction to "
            << featurePointStr().name() << endl;
    }


    // For points on multiple normals calculate intersection
    patchAttraction.setSize(pp.nPoints());
    patchAttraction = vector::zero;

    forAll(pp.localPoints(), pointI)
    {
        // Collect all different directions
        DynamicList<point> surfacePoints;
        DynamicList<vector> surfaceNormals;
        DynamicList<label> surfaceCount;

        binFeatureFaces
        (
            iter,
            featureCos,

            pp,
            snapDist,
            pointI,

            pointFaceNormals,
            pointFaceDisp,
            pointFaceCentres,

            surfacePoints,
            surfaceNormals,
            surfaceCount
        );

        const point& pt = pp.localPoints()[pointI];

        // Check the number of directions
        if (surfaceNormals.size() == 2)
        {
            plane pl0(surfacePoints[0], surfaceNormals[0]);
            plane pl1(surfacePoints[1], surfaceNormals[1]);
            plane::ray r(pl0.planeIntersect(pl1));
            vector n = r.dir() / mag(r.dir());

            // Get nearest point on infinite ray
            vector d = r.refPoint()-pt;
            d -= (d&n)*n;

            // Correct for attraction to non-dominant face
            correctAttraction
            (
                surfacePoints,
                surfaceCount,
                r.refPoint(),
                n,                  // normalised normal
                pt,

                d                   // perpendicular offset vector
            );

            // Trim to snap distance
            if (magSqr(d) > sqr(snapDist[pointI]))
            {
                d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
            }

            patchAttraction[pointI] = d;

            // Store constraints
            patchConstraints[pointI].applyConstraint(surfaceNormals[0]);
            patchConstraints[pointI].applyConstraint(surfaceNormals[1]);


            // Dump vector from point to points on faces
            if (featureEdgeStr.valid())
            {
                meshTools::writeOBJ(featureEdgeStr(), pt);
                featureEdgeVertI++;
                meshTools::writeOBJ(featureEdgeStr(), surfacePoints[0]);
                featureEdgeVertI++;
                meshTools::writeOBJ(featureEdgeStr(), surfacePoints[1]);
                featureEdgeVertI++;
                featureEdgeStr()
                    << "l " << featureEdgeVertI-2 << ' '
                    << featureEdgeVertI-1 << nl
                    << "l " << featureEdgeVertI-2 << ' '
                    << featureEdgeVertI << nl;
            }
        }
        else if (surfaceNormals.size() == 3)
        {
            // Calculate point from the faces.
            plane pl0(surfacePoints[0], surfaceNormals[0]);
            plane pl1(surfacePoints[1], surfaceNormals[1]);
            plane pl2(surfacePoints[2], surfaceNormals[2]);
            point cornerPt(pl0.planePlaneIntersect(pl1, pl2));

            vector d = cornerPt - pt;

            if (magSqr(d) > sqr(snapDist[pointI]))
            {
                d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
            }

            patchAttraction[pointI] = d;

            // Store constraints
            patchConstraints[pointI].applyConstraint(surfaceNormals[0]);
            patchConstraints[pointI].applyConstraint(surfaceNormals[1]);
            patchConstraints[pointI].applyConstraint(surfaceNormals[2]);

            // Dump vector from point to points on faces
            if (featurePointStr.valid())
            {
                meshTools::writeOBJ(featurePointStr(), pt);
                featurePointVertI++;
                meshTools::writeOBJ(featurePointStr(), surfacePoints[0]);
                featurePointVertI++;
                meshTools::writeOBJ(featurePointStr(), surfacePoints[1]);
                featurePointVertI++;
                meshTools::writeOBJ(featurePointStr(), surfacePoints[2]);
                featurePointVertI++;
                featurePointStr()
                    << "l " << featurePointVertI-3 << ' '
                    << featurePointVertI-2 << nl
                    << "l " << featurePointVertI-3 << ' '
                    << featurePointVertI-1 << nl
                    << "l " << featurePointVertI-3 << ' '
                    << featurePointVertI << nl;
            }
        }
    }

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        dumpMove
        (
            meshRefiner_.mesh().time().path()
          / "patchAttraction_" + name(iter) + ".obj",
            pp.localPoints(),
            pp.localPoints() + patchAttraction
        );
    }
}


// Finds nearest feature (within snapDist) for all points of pp.
void Foam::autoSnapDriver::determineAllFeatures
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const List<List<point> >& pointFaceNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,

    List<labelList>& pointAttractor,
    List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OFstream> featureEdgeStr;
    label featureEdgeVertI = 0;
    autoPtr<OFstream> featurePointStr;
    label featurePointVertI = 0;

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        featureEdgeStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featureEdge_" + name(iter) + ".obj"
            )
        );
        Pout<< nl << "Dumping feature-edge direction to "
            << featureEdgeStr().name() << endl;

        featurePointStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featurePoint_" + name(iter) + ".obj"
            )
        );
        Pout<< nl << "Dumping feature-point direction to "
            << featurePointStr().name() << endl;
    }

    const refinementFeatures& features = meshRefiner_.features();

    // Look at near feature edges
    labelList nearEdgeFeat;
    List<pointIndexHit> nearEdgeInfo;

    labelList nearPointFeat;
    labelList nearPointIndex;
    {
        scalarField snapDistSqr(sqr(snapDist));
        features.findNearestEdge
        (
            pp.localPoints(),
            snapDistSqr,
            nearEdgeFeat,
            nearEdgeInfo
        );

        // Look at near feature points
        features.findNearestPoint
        (
            pp.localPoints(),
            snapDistSqr,
            nearPointFeat,
            nearPointIndex
        );
    }

    forAll(pp.localPoints(), pointI)
    {
        const point& pt = pp.localPoints()[pointI];
        const label featI = nearEdgeFeat[pointI];
        const pointIndexHit& nearEdge = nearEdgeInfo[pointI];

        // Mark point on the nearest feature edge.
        if (nearEdge.hit())
        {
            // So we have a point on the feature edge. Use this instead
            // of our estimate from planes.

            label featEdgeI = nearEdge.index();
            edgeAttractors[featI][featEdgeI].append(nearEdge.hitPoint());
            const featureEdgeMesh& eMesh = features[nearEdgeFeat[pointI]];
            const edge& e = eMesh.edges()[featEdgeI];
            vector eVec = e.vec(eMesh.points());
            eVec /= mag(eVec);
            pointConstraint c;
            c.first() = 2;
            c.second() = eVec;
            edgeConstraints[featI][featEdgeI].append(c);

            // Store for later use
            patchAttraction[pointI] = nearEdge.hitPoint()-pt;
            patchConstraints[pointI] = c;

            // Dump
            if (featureEdgeStr.valid())
            {
                meshTools::writeOBJ(featureEdgeStr(), pt);
                featureEdgeVertI++;
                meshTools::writeOBJ(featureEdgeStr(), nearEdge.hitPoint());
                featureEdgeVertI++;
                featureEdgeStr()
                    << "l " << featureEdgeVertI-1 << ' '
                    << featureEdgeVertI << nl;
            }
        }


        // Mark point on the nearest feature point.
        if (nearPointFeat[pointI] != -1)
        {
            label featI = nearPointFeat[pointI];
            label index = nearPointIndex[pointI];
            const treeDataPoint& shapes = features.pointTrees()[featI].shapes();
            label featPointI = shapes.pointLabels()[index];
            const point& featPt = shapes.points()[featPointI];

            pointAttractor[featI][featPointI] = pointI;
            pointConstraints[featI][featPointI].first() = 3;
            pointConstraints[featI][featPointI].second() = vector::zero;

            // Dump
            if (featurePointStr.valid())
            {
                meshTools::writeOBJ(featurePointStr(), pt);
                featurePointVertI++;
                meshTools::writeOBJ(featurePointStr(), featPt);
                featurePointVertI++;
                featurePointStr()
                    << "l " << featurePointVertI-1 << ' '
                    << featurePointVertI << nl;
            }
        }
    }
}


void Foam::autoSnapDriver::determineFeatures
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    //const vectorField& faceSurfaceNormal,
    //const vectorField& faceDisp,
    //const vectorField& faceRotation,
    const List<List<point> >& pointFaceNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OFstream> featureEdgeStr;
    label featureEdgeVertI = 0;
    autoPtr<OFstream> missedEdgeStr;
    label missedVertI = 0;
    autoPtr<OFstream> featurePointStr;
    label featurePointVertI = 0;

    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        featureEdgeStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featureEdge_" + name(iter) + ".obj"
            )
        );
        Pout<< nl << "Dumping feature-edge sampling to "
            << featureEdgeStr().name() << endl;

        missedEdgeStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Pout<< nl << "Dumping feature-edges that are too far away to "
            << missedEdgeStr().name() << endl;

        featurePointStr.reset
        (
            new OFstream
            (
                meshRefiner_.mesh().time().path()
              / "featurePoint_" + name(iter) + ".obj"
            )
        );
        Pout<< nl << "Dumping feature-point sampling to "
            << featurePointStr().name() << endl;
    }

    const refinementFeatures& features = meshRefiner_.features();

    forAll(pp.localPoints(), pointI)
    {
        if (patchConstraints[pointI].first() == 0)
        {
            const point& pt = pp.localPoints()[pointI];

            // Collect all different directions
            DynamicList<point> surfacePoints;
            DynamicList<vector> surfaceNormals;
            DynamicList<label> surfaceCount;

            binFeatureFaces
            (
                iter,
                featureCos,

                pp,
                snapDist,
                pointI,

                //faceSurfaceNormal,
                //faceDisp,
                pointFaceNormals,
                pointFaceDisp,
                pointFaceCentres,

                surfacePoints,
                surfaceNormals,
                surfaceCount
            );

            // Check the number of directions
            if (surfaceNormals.size() == 2)
            {
                plane pl0(surfacePoints[0], surfaceNormals[0]);
                plane pl1(surfacePoints[1], surfaceNormals[1]);
                plane::ray r(pl0.planeIntersect(pl1));
                const vector n = r.dir() / mag(r.dir());

                // Get nearest point on infinite ray
                vector d = r.refPoint()-pt;
                d -= (d&n)*n;

                // Trim to snap distance
                if (magSqr(d) > sqr(snapDist[pointI]))
                {
                    d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
                }


                // Mark point on the nearest feature edge. Note that we
                // only search within the surrounding since the plane
                // reconstruction might find a feature where there isn't one.
                const point estimatedPt(pt + d);

                labelList nearEdgeFeat;
                List<pointIndexHit> nearEdgeInfo;
                features.findNearestEdge
                (
                    pointField(1, estimatedPt),
                    scalarField(1, sqr(snapDist[pointI])),
                    nearEdgeFeat,
                    nearEdgeInfo
                );

                const pointIndexHit& nearInfo = nearEdgeInfo[0];
                label featI = nearEdgeFeat[0];

                if (nearInfo.hit())
                {
                    // So we have a point on the feature edge. Use this instead
                    // of our estimate from planes.

                    edgeAttractors[featI][nearInfo.index()].append
                    (
                        nearInfo.hitPoint()
                    );
                    pointConstraint c;
                    c.applyConstraint(surfaceNormals[0]);
                    c.applyConstraint(surfaceNormals[1]);
                    edgeConstraints[featI][nearInfo.index()].append(c);

                    // Store for later use
                    patchAttraction[pointI] = nearInfo.hitPoint()-pt;
                    patchConstraints[pointI] = c;

                    // Dump
                    if (featureEdgeStr.valid())
                    {
                        meshTools::writeOBJ(featureEdgeStr(), pt);
                        featureEdgeVertI++;
                        meshTools::writeOBJ
                        (
                            featureEdgeStr(),
                            nearInfo.hitPoint()
                        );
                        featureEdgeVertI++;
                        featureEdgeStr()
                            << "l " << featureEdgeVertI-1 << ' '
                            << featureEdgeVertI << nl;
                    }
                }
                else
                {
                    if (missedEdgeStr.valid())
                    {
                        meshTools::writeOBJ(missedEdgeStr(), pt);
                        missedVertI++;
                        meshTools::writeOBJ
                        (
                            missedEdgeStr(),
                            nearInfo.missPoint()
                        );
                        missedVertI++;
                        missedEdgeStr()
                            << "l " << missedVertI-1 << ' '
                            << missedVertI << nl;
                    }
                }
            }
            else if (surfaceNormals.size() == 3)
            {
                // Calculate point from the faces.
                plane pl0(surfacePoints[0], surfaceNormals[0]);
                plane pl1(surfacePoints[1], surfaceNormals[1]);
                plane pl2(surfacePoints[2], surfaceNormals[2]);
                point cornerPt(pl0.planePlaneIntersect(pl1, pl2));
                vector d = cornerPt - pt;

                // Trim to snap distance
                if (magSqr(d) > sqr(snapDist[pointI]))
                {
                    d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
                }


                // Mark point on the nearest feature point.
                const point estimatedPt(pt + d);

                labelList nearPointFeat;
                labelList nearPointIndex;
                features.findNearestPoint
                (
                    pointField(1, estimatedPt),
                    scalarField(1, sqr(snapDist[pointI])),
                    nearPointFeat,
                    nearPointIndex
                );


                if (nearPointIndex[0] != -1)
                {
                    label featI = nearPointFeat[0];
                    label index = nearPointIndex[0];

                    const treeDataPoint& shapes =
                        features.pointTrees()[featI].shapes();
                    label featPointI = shapes.pointLabels()[index];
                    const point& featPt = shapes.points()[featPointI];

                    pointAttractor[featI][featPointI] = pointI;
                    pointConstraint& c = pointConstraints[featI][featPointI];
                    c.applyConstraint(surfaceNormals[0]);
                    c.applyConstraint(surfaceNormals[1]);
                    c.applyConstraint(surfaceNormals[2]);

                    // Store for later use
                    patchAttraction[pointI] = featPt-pt;
                    patchConstraints[pointI] = c;

                    // Dump
                    if (featurePointStr.valid())
                    {
                        meshTools::writeOBJ(featurePointStr(), pt);
                        featurePointVertI++;
                        meshTools::writeOBJ(featurePointStr(), featPt);
                        featurePointVertI++;
                        featurePointStr()
                            << "l " << featurePointVertI-1 << ' '
                            << featurePointVertI << nl;
                    }
                }
            }
        }
    }
}


void Foam::autoSnapDriver::featureAttractionUsingFeatureEdges
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    //const vectorField& faceSurfaceNormal,
    //const vectorField& faceDisp,
    //const vectorField& faceRotation,
    const List<List<point> >& pointFaceNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Collect ordered attractions on feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per feature, per feature-edge a list of attraction points and their
    // originating vertex.
    List<List<DynamicList<point> > > edgeAttractors(features.size());
    List<List<DynamicList<pointConstraint> > > edgeConstraints(features.size());
    forAll(features, featI)
    {
        label nFeatEdges = features[featI].edges().size();
        edgeAttractors[featI].setSize(nFeatEdges);
        edgeConstraints[featI].setSize(nFeatEdges);
    }

    // Per feature, per feature-point the pp point that is attracted to it.
    // This list is only used to subset the feature-points that are actually
    // used.
    List<labelList> pointAttractor(features.size());
    List<List<pointConstraint> > pointConstraints(features.size());
    forAll(features, featI)
    {
        label nFeatPoints = features[featI].points().size();
        pointAttractor[featI].setSize(nFeatPoints, -1);
        pointConstraints[featI].setSize(nFeatPoints);
    }

    // Reverse: from pp point to nearest feature
    vectorField allPatchAttraction(pp.nPoints(), vector::zero);
    List<pointConstraint> allPatchConstraints(pp.nPoints());

    determineFeatures
    //determineAllFeatures
    (
        iter,
        featureCos,

        pp,
        snapDist,

        pointFaceNormals,
        pointFaceDisp,
        pointFaceCentres,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,
        // pp point to nearest feature
        allPatchAttraction,
        allPatchConstraints
    );


    // Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Reverse lookup : go through all edgeAttractors and find the
    // nearest point on pp

    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint(pp.localPoints()),
        treeBoundBox(pp.localPoints()), // overall search domain
        8,                              // maxLevel
        10,                             // leafsize
        3.0                             // duplicity
    );

    // Per mesh point the point on nearest feature edge.
    patchAttraction.setSize(pp.nPoints());
    patchAttraction = vector::zero;
    patchConstraints.setSize(pp.nPoints());
    patchConstraints = pointConstraint();

    forAll(edgeAttractors, featI)
    {
        const List<DynamicList<point> >& edgeAttr = edgeAttractors[featI];
        const List<DynamicList<pointConstraint> >& edgeConstr =
            edgeConstraints[featI];

        forAll(edgeAttr, featEdgeI)
        {
            const DynamicList<point>& attr = edgeAttr[featEdgeI];
            forAll(attr, i)
            {
                // Find nearest pp point

                const point& featPt = attr[i];
                pointIndexHit nearInfo = ppTree.findNearest(featPt, sqr(GREAT));

                if (nearInfo.hit())
                {
                    label pointI = nearInfo.index();
                    const point attraction = featPt-pp.localPoints()[pointI];

                    // Check if this point is already being attracted. If so
                    // override it only if nearer.
                    if
                    (
                        patchConstraints[pointI].first() == 0
                     || magSqr(attraction) < magSqr(patchAttraction[pointI])
                    )
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = edgeConstr[featEdgeI][i];
                    }
                }
                else
                {
                    WarningIn
                    (
                        "autoSnapDriver::featureAttractionUsingFeatureEdges(..)"
                    )   << "Did not find pp point near " << featPt
                        << endl;
                }
            }
        }
    }


    // Different procs might have different patchAttraction,patchConstraints
    // however these only contain geometric information, no topology
    // so as long as we synchronise after overriding with feature points
    // there is no problem, just possibly a small error.


    // Find nearest mesh point to feature point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (overrides attraction to feature edge)
    forAll(pointAttractor, featI)
    {
        const labelList& pointAttr = pointAttractor[featI];
        const List<pointConstraint>& pointConstr = pointConstraints[featI];

        forAll(pointAttr, featPointI)
        {
            if (pointAttr[featPointI] != -1)
            {
                const point& featPt = features[featI].points()
                [
                    featPointI
                ];

                // Find nearest pp point
                pointIndexHit nearInfo = ppTree.findNearest(featPt, sqr(GREAT));

                if (nearInfo.hit())
                {
                    label pointI = nearInfo.index();
                    const point& pt = pp.localPoints()[pointI];
                    const point attraction = featPt-pt;

                    // - already attracted to feature edge : point always wins
                    // - already attracted to feature point: nearest wins

                    if (patchConstraints[pointI].first() <= 1)
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = pointConstr[featPointI];
                    }
                    else if (patchConstraints[pointI].first() == 2)
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = pointConstr[featPointI];
                    }
                    else if (patchConstraints[pointI].first() == 3)
                    {
                        // Only if nearer
                        if
                        (
                            magSqr(attraction)
                          < magSqr(patchAttraction[pointI])
                        )
                        {
                            patchAttraction[pointI] = attraction;
                            patchConstraints[pointI] = pointConstr[featPointI];
                        }
                    }
                }
            }
        }
    }


    // Dump
    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        OFstream featureEdgeStr
        (
            meshRefiner_.mesh().time().path()
          / "edgeAttractors_" + name(iter) + ".obj"
        );
        label featureEdgeVertI = 0;
        Pout<< nl << "Dumping feature-edge attraction to "
            << featureEdgeStr.name() << endl;

        OFstream featurePointStr
        (
            meshRefiner_.mesh().time().path()
          / "pointAttractors_" + name(iter) + ".obj"
        );
        label featurePointVertI = 0;
        Pout<< nl << "Dumping feature-point attraction to "
            << featurePointStr.name() << endl;

        forAll(patchConstraints, pointI)
        {
            const point& pt = pp.localPoints()[pointI];

            if (patchConstraints[pointI].first() == 2)
            {
                meshTools::writeOBJ(featureEdgeStr, pt);
                featureEdgeVertI++;
                meshTools::writeOBJ(featureEdgeStr, pt+patchAttraction[pointI]);
                featureEdgeVertI++;
                featureEdgeStr << "l " << featureEdgeVertI-1
                    << ' ' << featureEdgeVertI << nl;
            }
            else if (patchConstraints[pointI].first() == 3)
            {
                meshTools::writeOBJ(featurePointStr, pt);
                featurePointVertI++;
                meshTools::writeOBJ
                (
                    featurePointStr,
                    pt+patchAttraction[pointI]
                );
                featurePointVertI++;
                featurePointStr << "l " << featurePointVertI-1
                    << ' ' << featurePointVertI << nl;
            }
        }
    }



    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in allPatchConstraints)

    while (true)
    {
        label nChanged = 0;

        const labelListList& pointEdges = pp.pointEdges();
        forAll(pointEdges, pointI)
        {
            if (patchConstraints[pointI].first() == 2)
            {
                const point& pt = pp.localPoints()[pointI];
                const labelList& pEdges = pointEdges[pointI];
                const vector& featVec = patchConstraints[pointI].second();

                // Detect whether there are edges in both directions.
                // (direction along the feature edge that is)
                bool hasPos = false;
                bool hasNeg = false;

                forAll(pEdges, pEdgeI)
                {
                    const edge& e = pp.edges()[pEdges[pEdgeI]];
                    label nbrPointI = e.otherVertex(pointI);

                    if (patchConstraints[nbrPointI].first() != 0)
                    {
                        const point& nbrPt = pp.localPoints()[nbrPointI];
                        const point featPt = nbrPt + patchAttraction[nbrPointI];
                        const scalar cosAngle = (featVec & (featPt-pt));

                        if (cosAngle > 0)
                        {
                            hasPos = true;
                        }
                        else
                        {
                            hasNeg = true;
                        }
                    }
                }

                if (!hasPos || !hasNeg)
                {
                    //Pout<< "**Detected feature string end at  "
                    //    << pp.localPoints()[pointI] << endl;

                    // No string. Assign best choice on either side
                    label bestPosPointI = -1;
                    scalar minPosDistSqr = GREAT;
                    label bestNegPointI = -1;
                    scalar minNegDistSqr = GREAT;

                    forAll(pEdges, pEdgeI)
                    {
                        const edge& e = pp.edges()[pEdges[pEdgeI]];
                        label nbrPointI = e.otherVertex(pointI);

                        if
                        (
                            patchConstraints[nbrPointI].first() == 0
                         && allPatchConstraints[nbrPointI].first() != 0
                        )
                        {
                            const vector& nbrFeatVec =
                                allPatchConstraints[pointI].second();

                            if (mag(featVec&nbrFeatVec) > featureCos)
                            {
                                // nbrPointI attracted to sameish feature
                                // Note: also check on position.

                                scalar d2 = magSqr
                                (
                                    allPatchAttraction[nbrPointI]
                                );

                                const point featPt =
                                    pp.localPoints()[nbrPointI]
                                  + allPatchAttraction[nbrPointI];
                                const scalar cosAngle = (featVec & (featPt-pt));

                                if (cosAngle > 0)
                                {
                                    if (!hasPos && d2 < minPosDistSqr)
                                    {
                                        minPosDistSqr = d2;
                                        bestPosPointI = nbrPointI;
                                    }
                                }
                                else
                                {
                                    if (!hasNeg && d2 < minNegDistSqr)
                                    {
                                        minNegDistSqr = d2;
                                        bestNegPointI = nbrPointI;
                                    }
                                }
                            }
                        }
                    }

                    if (bestPosPointI != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt = pp.localPoints()[bestPosPointI];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << allPatchAttraction[bestPosPointI]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestPosPointI] =
                            0.5*allPatchAttraction[bestPosPointI];
                        patchConstraints[bestPosPointI] =
                            allPatchConstraints[bestPosPointI];

                        nChanged++;
                    }
                    if (bestNegPointI != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt = pp.localPoints()[bestNegPointI];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << allPatchAttraction[bestNegPointI]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestNegPointI] =
                            0.5*allPatchAttraction[bestNegPointI];
                        patchConstraints[bestNegPointI] =
                            allPatchConstraints[bestNegPointI];

                        nChanged++;
                    }
                }
            }
        }


        reduce(nChanged, sumOp<label>());
        Info<< "Stringing feature edges : changed " << nChanged << " points"
            << endl;
        if (nChanged == 0)
        {
            break;
        }
    }



    // Avoid diagonal attraction
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Attract one of the non-diagonal points.

    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];
        // For now just detect any attraction. Improve this to look at
        // actual attraction position and only if would form a problem add
        // the non-diagonal point
        if (f.size() == 4)
        {
            label nAttract = 0;
            label firstAttract = -1;
            forAll(f, fp)
            {
                label pointI = f[fp];
                if (patchConstraints[pointI].first() == 2)
                {
                    nAttract++;
                    if (firstAttract == -1)
                    {
                        firstAttract = fp;
                    }
                }
            }
            if (nAttract == 2)
            {
                label nextAttract = f.fcIndex(f.fcIndex(firstAttract));
                label pointI = f[nextAttract];

                if (patchConstraints[pointI].first() == 2)
                {
                    // Found two diagonal points that being attracted.
                    // For now just attract my one to the average of those.
                    const label i0 = f[firstAttract];
                    const point pt0 = pp.localPoints()[i0]+patchAttraction[i0];
                    const label i1 = f[nextAttract];
                    const point pt1 = pp.localPoints()[i1]+patchAttraction[i1];
                    const point mid = 0.5*(pt0+pt1);


                    const scalar cosAngle = mag
                    (
                        patchConstraints[i0].second()
                      & patchConstraints[i1].second()
                    );

                    //Pout<< "Found diagonal attraction at indices:"
                    //    << firstAttract
                    //    << " and " << nextAttract
                    //    << " with cosAngle:" << cosAngle
                    //    << " mid:" << mid << endl;

                    if (cosAngle > featureCos)
                    {
                        // Add the nearest of the other two points as attractor
                        label minFp = -1;
                        scalar minDistSqr = GREAT;
                        forAll(f, fp)
                        {
                            label pointI = f[fp];
                            if (patchConstraints[pointI].first() == 0)
                            {
                                const point& pt = pp.localPoints()[pointI];
                                scalar distSqr = magSqr(mid-pt);
                                if (distSqr < minDistSqr)
                                {
                                    distSqr = minDistSqr;
                                    minFp = fp;
                                }
                            }
                        }
                        if (minFp != -1)
                        {
                            label minPointI = f[minFp];
                            patchAttraction[minPointI] =
                                mid-pp.localPoints()[minPointI];
                            patchConstraints[minPointI] =
                                patchConstraints[f[firstAttract]];
                        }
                    }
                    else
                    {
                        //Pout<< "Diagonal attractors at" << nl
                        //    << "    pt0:" << pt0
                        //    << "    constraint:"
                        //    << patchConstraints[i0].second() << nl
                        //    << "    pt1:" << pt1
                        //    << "    constraint:"
                        //    << patchConstraints[i1].second() << nl
                        //    << "    make too large an angle:"
                        //    <<  mag
                        //        (
                        //            patchConstraints[i0].second()
                        //          & patchConstraints[i1].second()
                        //        )
                        //    << endl;
                    }
                }
            }
        }
    }


    if (debug&meshRefinement::OBJINTERSECTIONS)
    {
        dumpMove
        (
            meshRefiner_.mesh().time().path()
          / "patchAttraction_" + name(iter) + ".obj",
            pp.localPoints(),
            pp.localPoints() + patchAttraction
        );
    }
}


// Correct for squeezing of face
void Foam::autoSnapDriver::preventFaceSqueeze
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    pointField points;
    face singleF;
    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        if (f.size() != points.size())
        {
            points.setSize(f.size());
            singleF.setSize(f.size());
            for (label i = 0; i < f.size(); i++)
            {
                singleF[i] = i;
            }
        }
        label nConstraints = 0;
        forAll(f, fp)
        {
            label pointI = f[fp];
            if (patchConstraints[pointI].first() != 0)
            {
                points[fp] = pp.localPoints()[pointI] + patchAttraction[pointI];
                nConstraints++;
            }
            else
            {
                points[fp] = pp.localPoints()[pointI];
            }
        }

        if (nConstraints == f.size())
        {
            scalar oldArea = f.mag(pp.localPoints());
            scalar newArea = singleF.mag(points);
            if (newArea < 0.1*oldArea)
            {
                // For now remove the point with largest distance
                label maxFp = -1;
                scalar maxS = -1;
                forAll(f, fp)
                {
                    scalar s = magSqr(patchAttraction[f[fp]]);
                    if (s > maxS)
                    {
                        maxS = s;
                        maxFp = fp;
                    }
                }
                if (maxFp != -1)
                {
                    label pointI = f[maxFp];
                    // Lower attraction on pointI
                    patchAttraction[pointI] *= 0.5;
                }
            }
        }
    }
}


Foam::vectorField Foam::autoSnapDriver::calcNearestSurfaceFeature
(
    const label iter,
    const scalar featureCos,
    const scalar featureAttract,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    motionSmoother& meshMover
) const
{
    Info<< "Calculating patchDisplacement as distance to nearest surface"
        << " point ..." << endl;

    const indirectPrimitivePatch& pp = meshMover.patch();
    const pointField& localPoints = pp.localPoints();
    const fvMesh& mesh = meshRefiner_.mesh();

    // Displacement and orientation per pp face.
    vectorField faceDisp(pp.size(), vector::zero);
    vectorField faceSurfaceNormal(pp.size(), vector::zero);
    vectorField faceRotation(pp.size(), vector::zero);

    calcNearestFace
    (
        iter,
        pp,
        faceDisp,
        faceSurfaceNormal,
        faceRotation
    );

    // Start off with nearest point on surface.
    vectorField patchDisp = nearestDisp;


    // Collect (possibly remote) face-wise data on coupled points.
    // - faceSurfaceNormal
    // - faceDisp
    // - faceRotation
    // - faceCentres

    // For now just get all surrounding face data. Expensive - should just
    // store and sync data on coupled points only (see e.g PatchToolsNormals.C)

    List<List<point> > pointFaceNormals(pp.nPoints());
    List<List<point> > pointFaceDisp(pp.nPoints());
    List<List<point> > pointFaceCentres(pp.nPoints());

    // Fill local data
    forAll(pp.pointFaces(), pointI)
    {
        const labelList& pFaces = pp.pointFaces()[pointI];
        List<point>& pNormals = pointFaceNormals[pointI];
        pNormals.setSize(pFaces.size());
        List<point>& pDisp = pointFaceDisp[pointI];
        pDisp.setSize(pFaces.size());
        List<point>& pFc = pointFaceCentres[pointI];
        pFc.setSize(pFaces.size());
        forAll(pFaces, i)
        {
            pNormals[i] = faceSurfaceNormal[pFaces[i]];
            pDisp[i] = faceDisp[pFaces[i]];
            pFc[i] = pp.faceCentres()[pFaces[i]];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceNormals,
        listPlusEqOp(),
        List<point>(),
        listTransform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceDisp,
        listPlusEqOp(),
        List<point>(),
        listTransform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceCentres,
        listPlusEqOp(),
        List<point>(),
        listTransform()
    );



    // Nearest feature
    vectorField patchAttraction(localPoints.size(), vector::zero);
    // Constraints at feature
    List<pointConstraint> patchConstraints(localPoints.size());
    //featureAttractionUsingReconstruction
    featureAttractionUsingFeatureEdges
    (
        iter,
        featureCos,

        pp,
        snapDist,

        pointFaceNormals,
        pointFaceDisp,
        pointFaceCentres,

        patchAttraction,
        patchConstraints
    );

    preventFaceSqueeze
    (
        iter,
        featureCos,

        pp,
        snapDist,

        patchAttraction,
        patchConstraints
    );


    Info<< "Attraction:" << endl
        << "     linear   : max:" << gMax(patchDisp)
        << " avg:" << gAverage(patchDisp)
        << endl
        << "     feature  : max:" << gMax(patchAttraction)
        << " avg:" << gAverage(patchAttraction)
        << endl;


    // So now we have:
    // - patchDisp          : point movement to go to nearest point on surface
    //                       (either direct or through interpolation of
    //                        face nearest)
    // - patchAttraction    : direct attraction to features
    // - patchConstraints   : type of features

    // Use any combination of patchDisp and direct feature
    // attraction.


    // Mix with direct feature attraction
    forAll(patchConstraints, pointI)
    {
        if (patchConstraints[pointI].first() != 0)
        {
            patchDisp[pointI] =
                (1.0-featureAttract)*patchDisp[pointI]
              + featureAttract*patchAttraction[pointI];
        }
    }

    //dumpMove
    //(
    //    mesh.time().path()
    //  / "linearPatchDisp_" + name(iter) + ".obj",
    //    pp.localPoints(),
    //    pp.localPoints() + patchDisp
    //);



    // Count
    {
        label nPlanar = 0;
        label nEdge = 0;
        label nPoint = 0;

        forAll(patchConstraints, pointI)
        {
            if (patchConstraints[pointI].first() == 1)
            {
                nPlanar++;
            }
            else if (patchConstraints[pointI].first() == 2)
            {
                nEdge++;
            }
            else if (patchConstraints[pointI].first() == 3)
            {
                nPoint++;
            }
        }

        label nTotPoints = returnReduce(pp.nPoints(), sumOp<label>());
        reduce(nPlanar, sumOp<label>());
        reduce(nEdge, sumOp<label>());
        reduce(nPoint, sumOp<label>());
        Info<< "Feature analysis : total points:"
            << nTotPoints
            << " attraction to :" << nl
            << "    feature point   : " << nPoint << nl
            << "    feature edge    : " << nEdge << nl
            << "    nearest surface : " << nTotPoints-nPoint-nEdge
            << " (rest)" << nl
            << endl;
    }



    // Now we have the displacement per patch point to move onto the surface
    // Split into tangential and normal direction.
    // - start off with all non-constrained points following the constrained
    //   ones since point normals not relevant.
    // - finish with only tangential component smoothed.
    // Note: tangential is most
    // likely to come purely from face-centre snapping, not face rotation.
    if (featureAttract < 1-0.001)
    {
        // 1. Smoothed all displacement
        vectorField smoothedPatchDisp = patchDisp;
        smoothAndConstrain
        (
            pp,
            patchConstraints,
            smoothedPatchDisp
        );


        // 2. Smoothed tangential component
        vectorField tangPatchDisp = patchDisp;
        tangPatchDisp -= (pp.pointNormals() & patchDisp) * pp.pointNormals();
        smoothAndConstrain
        (
            pp,
            patchConstraints,
            tangPatchDisp
        );

        // Re-add normal component
        tangPatchDisp += (pp.pointNormals() & patchDisp) * pp.pointNormals();

        if (debug&meshRefinement::OBJINTERSECTIONS)
        {
            dumpMove
            (
                mesh.time().path()
              / "tangPatchDispConstrained_" + name(iter) + ".obj",
                pp.localPoints(),
                pp.localPoints() + tangPatchDisp
            );
        }

        patchDisp =
             (1.0-featureAttract)*smoothedPatchDisp
           + featureAttract*tangPatchDisp;
    }



    const scalar relax = featureAttract; //1.0;
    patchDisp *= relax;


    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value (note: cannot use VGREAT)
    );

    return patchDisp;
}


// ************************************************************************* //
