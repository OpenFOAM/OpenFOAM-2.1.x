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

Application
    surfaceFeatureExtract

Description
    Extracts and writes surface features to file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurface.H"
#include "surfaceFeatures.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"
#include "triSurfaceFields.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "unitConversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dumpBox(const treeBoundBox& bb, const fileName& fName)
{
    OFstream str(fName);

    Info<< "Dumping bounding box " << bb << " as lines to obj file "
        << str.name() << endl;


    pointField boxPoints(bb.points());

    forAll(boxPoints, i)
    {
        meshTools::writeOBJ(str, boxPoints[i]);
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str<< "l " << e[0]+1 <<  ' ' << e[1]+1 << nl;
    }
}


// Deletes all edges inside/outside bounding box from set.
void deleteBox
(
    const triSurface& surf,
    const treeBoundBox& bb,
    const bool removeInside,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    forAll(edgeStat, edgeI)
    {
        const point eMid = surf.edges()[edgeI].centre(surf.localPoints());

        if (removeInside ? bb.contains(eMid) : !bb.contains(eMid))
        {
            edgeStat[edgeI] = surfaceFeatures::NONE;
        }
    }
}


void drawHitProblem
(
    label fI,
    const triSurface& surf,
    const pointField& start,
    const pointField& faceCentres,
    const pointField& end,
    const List<pointIndexHit>& hitInfo
)
{
    Info<< nl << "# findLineAll did not hit its own face."
        << nl << "# fI " << fI
        << nl << "# start " << start[fI]
        << nl << "# f centre " << faceCentres[fI]
        << nl << "# end " << end[fI]
        << nl << "# hitInfo " << hitInfo
        << endl;

    meshTools::writeOBJ(Info, start[fI]);
    meshTools::writeOBJ(Info, faceCentres[fI]);
    meshTools::writeOBJ(Info, end[fI]);

    Info<< "l 1 2 3" << endl;

    meshTools::writeOBJ(Info, surf.points()[surf[fI][0]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][1]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][2]]);

    Info<< "f 4 5 6" << endl;

    forAll(hitInfo, hI)
    {
        label hFI = hitInfo[hI].index();

        meshTools::writeOBJ(Info, surf.points()[surf[hFI][0]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][1]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][2]]);

        Info<< "f "
            << 3*hI + 7 << " "
            << 3*hI + 8 << " "
            << 3*hI + 9
            << endl;
    }
}


// Unmark non-manifold edges if individual triangles are not features
void unmarkBaffles
(
    const triSurface& surf,
    const scalar includedAngle,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

    const labelListList& edgeFaces = surf.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            label i0 = eFaces[0];
            //const labelledTri& f0 = surf[i0];
            const Foam::vector& n0 = surf.faceNormals()[i0];

            //Pout<< "edge:" << edgeI << " n0:" << n0 << endl;

            bool same = true;

            for (label i = 1; i < eFaces.size(); i++)
            {
                //const labelledTri& f = surf[i];
                const Foam::vector& n = surf.faceNormals()[eFaces[i]];

                //Pout<< "    mag(n&n0): " << mag(n&n0) << endl;

                if (mag(n&n0) < minCos)
                {
                    same = false;
                    break;
                }
            }

            if (same)
            {
                edgeStat[edgeI] = surfaceFeatures::NONE;
            }
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract and write surface features to file"
    );
    argList::noParallel();
    argList::validArgs.append("surface");
    argList::validArgs.append("output set");

    argList::addOption
    (
        "includedAngle",
        "degrees",
        "construct feature set from included angle [0..180]"
    );
    argList::addOption
    (
        "set",
        "name",
        "use existing feature set from file"
    );
    argList::addOption
    (
        "minLen",
        "scalar",
        "remove features shorter than the specified cumulative length"
    );
    argList::addOption
    (
        "minElem",
        "int",
        "remove features with fewer than the specified number of edges"
    );
    argList::addOption
    (
        "subsetBox",
        "((x0 y0 z0)(x1 y1 z1))",
        "remove edges outside specified bounding box"
    );
    argList::addOption
    (
        "deleteBox",
        "((x0 y0 z0)(x1 y1 z1))",
        "remove edges within specified bounding box"
    );
    argList::addBoolOption
    (
        "writeObj",
        "write extendedFeatureEdgeMesh obj files"
    );
    argList::addBoolOption
    (
        "writeVTK",
        "write extendedFeatureEdgeMesh vtk files"
    );
    argList::addOption
    (
        "closeness",
        "scalar",
        "span to look for surface closeness"
    );
    argList::addOption
    (
        "featureProximity",
        "scalar",
        "distance to look for close features"
    );
    argList::addBoolOption
    (
        "manifoldEdgesOnly",
        "remove any non-manifold (open or more than two connected faces) edges"
    );

#   ifdef ENABLE_CURVATURE
    argList::addBoolOption
    (
        "calcCurvature",
        "calculate curvature and closeness fields"
    );
#   endif


#   include "setRootCase.H"
#   include "createTime.H"

    bool writeVTK = args.optionFound("writeVTK");

    bool writeObj = args.optionFound("writeObj");

    bool curvature = args.optionFound("curvature");

    if (curvature && env("FOAM_SIGFPE"))
    {
        WarningIn(args.executable())
            << "Detected floating point exception trapping (FOAM_SIGFPE)."
            << " This might give" << nl
            << "    problems when calculating curvature on straight angles"
            << " (infinite curvature)" << nl
            << "    Switch it off in case of problems." << endl;
    }


    Info<< "Feature line extraction is only valid on closed manifold surfaces."
        << endl;

    const fileName surfFileName = args[1];
    const fileName outFileName  = args[2];

    Info<< "Surface            : " << surfFileName << nl
        << "Output feature set : " << outFileName << nl
        << endl;

    fileName sFeatFileName = surfFileName.lessExt().name();


    // Read
    // ~~~~

    triSurface surf(surfFileName);

    Info<< "Statistics:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    faceList faces(surf.size());

    forAll(surf, fI)
    {
        faces[fI] = surf[fI].triFaceFace();
    }

    // Either construct features from surface&featureangle or read set.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    surfaceFeatures set(surf);

    if (args.optionFound("set"))
    {
        const fileName setName = args["set"];

        Info<< "Reading existing feature set from file " << setName << endl;

        set = surfaceFeatures(surf, setName);
    }
    else if (args.optionFound("includedAngle"))
    {
        const scalar includedAngle = args.optionRead<scalar>("includedAngle");

        Info<< "Constructing feature set from included angle " << includedAngle
            << endl;

        set = surfaceFeatures(surf, includedAngle);

        // Info<< nl << "Writing initial features" << endl;
        // set.write("initial.fSet");
        // set.writeObj("initial");
    }
    else
    {
        FatalErrorIn(args.executable())
            << "No initial feature set. Provide either one"
            << " of -set (to read existing set)" << nl
            << " or -includedAngle (to new set construct from angle)"
            << exit(FatalError);
    }

    Info<< nl
        << "Initial feature set:" << nl
        << "    feature points : " << set.featurePoints().size() << nl
        << "    feature edges  : " << set.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << set.nRegionEdges() << nl
        << "        external edges : " << set.nExternalEdges() << nl
        << "        internal edges : " << set.nInternalEdges() << nl
        << endl;


    // Trim set
    // ~~~~~~~~

    scalar minLen = -GREAT;
    if (args.optionReadIfPresent("minLen", minLen))
    {
        Info<< "Removing features of length < " << minLen << endl;
    }

    label minElem = 0;
    if (args.optionReadIfPresent("minElem", minElem))
    {
        Info<< "Removing features with number of edges < " << minElem << endl;
    }

    // Trim away small groups of features
    if (minElem > 0 || minLen > 0)
    {
        set.trimFeatures(minLen, minElem);
        Info<< endl << "Removed small features" << endl;
    }


    // Subset
    // ~~~~~~

    // Convert to marked edges, points
    List<surfaceFeatures::edgeStatus> edgeStat(set.toStatus());

    if (args.optionFound("subsetBox"))
    {
        treeBoundBox bb
        (
            args.optionLookup("subsetBox")()
        );

        Info<< "Removing all edges outside bb " << bb << endl;
        dumpBox(bb, "subsetBox.obj");

        deleteBox
        (
            surf,
            bb,
            false,
            edgeStat
        );
    }
    else if (args.optionFound("deleteBox"))
    {
        treeBoundBox bb
        (
            args.optionLookup("deleteBox")()
        );

        Info<< "Removing all edges inside bb " << bb << endl;
        dumpBox(bb, "deleteBox.obj");

        deleteBox
        (
            surf,
            bb,
            true,
            edgeStat
        );
    }

    if (args.optionFound("manifoldEdgesOnly"))
    {
        Info<< "Removing all non-manifold edges" << endl;

        forAll(edgeStat, edgeI)
        {
            if (surf.edgeFaces()[edgeI].size() != 2)
            {
                edgeStat[edgeI] = surfaceFeatures::NONE;
            }
        }
    }


    surfaceFeatures newSet(surf);
    newSet.setFromStatus(edgeStat);

    //Info<< endl << "Writing trimmed features to " << outFileName << endl;
    //newSet.write(outFileName);

    // Info<< endl << "Writing edge objs." << endl;
    // newSet.writeObj("final");

    Info<< nl
        << "Final feature set:" << nl
        << "    feature points : " << newSet.featurePoints().size() << nl
        << "    feature edges  : " << newSet.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << newSet.nRegionEdges() << nl
        << "        external edges : " << newSet.nExternalEdges() << nl
        << "        internal edges : " << newSet.nInternalEdges() << nl
        << endl;

    // Extracting and writing a extendedFeatureEdgeMesh

    extendedFeatureEdgeMesh feMesh
    (
        newSet,
        runTime,
        sFeatFileName + ".extendedFeatureEdgeMesh"
    );

    Info<< nl << "Writing extendedFeatureEdgeMesh to " << feMesh.objectPath()
        << endl;

    if (writeObj)
    {
        feMesh.writeObj(surfFileName.lessExt().name());
    }

    feMesh.write();


    // Write a featureEdgeMesh for backwards compatibility
    {
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                surfFileName.lessExt().name() + ".eMesh",   // name
                runTime.constant(), // instance
                "triSurface",
                runTime,            // registry
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.objectPath() << endl;

        bfeMesh.regIOobject::write();
    }

    triSurfaceMesh searchSurf
    (
        IOobject
        (
            sFeatFileName + ".closeness",
            runTime.constant(),
            "extendedFeatureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf
    );

    if (!curvature)
    {
        Info<< "End\n" << endl;

        return 0;
    }

    // Find close features

    // // Dummy trim operation to mark features
    // labelList featureEdgeIndexing = newSet.trimFeatures(-GREAT, 0);

    // scalarField surfacePtFeatureIndex(surf.points().size(), -1);

    // forAll(newSet.featureEdges(), eI)
    // {
    //     const edge& e = surf.edges()[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.start()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.end()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];
    // }

    // if (writeVTK)
    // {
    //     vtkSurfaceWriter().write
    //     (
    //         runTime.constant()/"triSurface",    // outputDir
    //         sFeatFileName,                      // surfaceName
    //         surf.points(),
    //         faces,
    //         "surfacePtFeatureIndex",            // fieldName
    //         surfacePtFeatureIndex,
    //         true,                               // isNodeValues
    //         true                                // verbose
    //     );
    // }

    // Random rndGen(343267);

    // treeBoundBox surfBB
    // (
    //     treeBoundBox(searchSurf.bounds()).extend(rndGen, 1e-4)
    // );

    // surfBB.min() -= Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    // surfBB.max() += Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // indexedOctree<treeDataEdge> ftEdTree
    // (
    //     treeDataEdge
    //     (
    //         false,
    //         surf.edges(),
    //         surf.localPoints(),
    //         newSet.featureEdges()
    //     ),
    //     surfBB,
    //     8,      // maxLevel
    //     10,     // leafsize
    //     3.0     // duplicity
    // );

    // labelList nearPoints = ftEdTree.findBox
    // (
    //     treeBoundBox
    //     (
    //         sPt - featureSearchSpan*Foam::vector::one,
    //         sPt + featureSearchSpan*Foam::vector::one
    //     )
    // );

    Info<< "Examine curvature, feature proximity and internal and "
        << "external closeness." << endl;

    // Internal and external closeness

    // Prepare start and end points for intersection tests

    const vectorField& normals = searchSurf.faceNormals();

    scalar span = searchSurf.bounds().mag();

    args.optionReadIfPresent("closeness", span);

    scalar externalAngleTolerance = 10;
    scalar externalToleranceCosAngle = Foam::cos
    (
        degToRad(180 - externalAngleTolerance)
    );

    scalar internalAngleTolerance = 45;
    scalar internalToleranceCosAngle = Foam::cos
    (
        degToRad(180 - internalAngleTolerance)
    );

    Info<< "externalToleranceCosAngle: " << externalToleranceCosAngle << nl
        << "internalToleranceCosAngle: " << internalToleranceCosAngle
        << endl;

    // Info<< "span " << span << endl;

    pointField start = searchSurf.faceCentres() - span*normals;
    pointField end = searchSurf.faceCentres() + span*normals;
    const pointField& faceCentres = searchSurf.faceCentres();

    List<List<pointIndexHit> > allHitInfo;

    // Find all intersections (in order)
    searchSurf.findLineAll(start, end, allHitInfo);

    scalarField internalCloseness(start.size(), GREAT);
    scalarField externalCloseness(start.size(), GREAT);

    forAll(allHitInfo, fI)
    {
        const List<pointIndexHit>& hitInfo = allHitInfo[fI];

        if (hitInfo.size() < 1)
        {
            drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

            // FatalErrorIn(args.executable())
            //     << "findLineAll did not hit its own face."
            //     << exit(FatalError);
        }
        else if (hitInfo.size() == 1)
        {
            if (!hitInfo[0].hit())
            {
                // FatalErrorIn(args.executable())
                //     << "findLineAll did not hit any face."
                //     << exit(FatalError);
            }
            else if (hitInfo[0].index() != fI)
            {
                drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

                // FatalErrorIn(args.executable())
                //     << "findLineAll did not hit its own face."
                //     << exit(FatalError);
            }
        }
        else
        {
            label ownHitI = -1;

            forAll(hitInfo, hI)
            {
                // Find the hit on the triangle that launched the ray

                if (hitInfo[hI].index() == fI)
                {
                    ownHitI = hI;

                    break;
                }
            }

            if (ownHitI < 0)
            {
                drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

                // FatalErrorIn(args.executable())
                //     << "findLineAll did not hit its own face."
                //     << exit(FatalError);
            }
            else if (ownHitI == 0)
            {
                // There are no internal hits, the first hit is the closest
                // external hit

                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI + 1].index()])
                  < externalToleranceCosAngle
                )
                {
                    externalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI + 1].hitPoint()
                    );
                }
            }
            else if (ownHitI == hitInfo.size() - 1)
            {
                // There are no external hits, the last but one hit is the
                // closest internal hit

                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI - 1].index()])
                  < internalToleranceCosAngle
                )
                {
                    internalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI - 1].hitPoint()
                    );
                }
            }
            else
            {
                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI + 1].index()])
                  < externalToleranceCosAngle
                )
                {
                    externalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI + 1].hitPoint()
                    );
                }

                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI - 1].index()])
                  < internalToleranceCosAngle
                )
                {
                    internalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI - 1].hitPoint()
                    );
                }
            }
        }
    }

    triSurfaceScalarField internalClosenessField
    (
        IOobject
        (
            sFeatFileName + ".internalCloseness",
            runTime.constant(),
            "extendedFeatureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        internalCloseness
    );

    internalClosenessField.write();

    triSurfaceScalarField externalClosenessField
    (
        IOobject
        (
            sFeatFileName + ".externalCloseness",
            runTime.constant(),
            "extendedFeatureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        externalCloseness
    );

    externalClosenessField.write();


#ifdef ENABLE_CURVATURE
    scalarField k = calcCurvature(surf);

    // Modify the curvature values on feature edges and points to be zero.

    forAll(newSet.featureEdges(), fEI)
    {
        const edge& e = surf.edges()[newSet.featureEdges()[fEI]];

        k[surf.meshPoints()[e.start()]] = 0.0;
        k[surf.meshPoints()[e.end()]] = 0.0;
    }

    triSurfacePointScalarField kField
    (
        IOobject
        (
            sFeatFileName + ".curvature",
            runTime.constant(),
            "extendedFeatureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        k
    );

    kField.write();
#endif

    if (writeVTK)
    {
        vtkSurfaceWriter().write
        (
            runTime.constant()/"triSurface",    // outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "internalCloseness",                // fieldName
            internalCloseness,
            false,                              // isNodeValues
            true                                // verbose
        );

        vtkSurfaceWriter().write
        (
            runTime.constant()/"triSurface",    // outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "externalCloseness",                // fieldName
            externalCloseness,
            false,                              // isNodeValues
            true                                // verbose
        );

#       ifdef ENABLE_CURVATURE
        vtkSurfaceWriter().write
        (
            runTime.constant()/"triSurface",    // outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "curvature",                        // fieldName
            k,
            true,                               // isNodeValues
            true                                // verbose
        );
#       endif
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
