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

Application
    snappyHexMesh

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "autoRefineDriver.H"
#include "autoSnapDriver.H"
#include "autoLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Check writing tolerance before doing any serious work
scalar getMergeDistance(const polyMesh& mesh, const scalar mergeTol)
{
    const boundBox& meshBb = mesh.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstream::ASCII)
    {
        const scalar writeTol = std::pow
        (
            scalar(10.0),
            -scalar(IOstream::defaultPrecision())
        );

        if (mergeTol < writeTol)
        {
            FatalErrorIn("getMergeDistance(const polyMesh&, const dictionary&)")
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


// Write mesh and additional information
void writeMesh
(
    const string& msg,
    const meshRefinement& meshRefiner,
    const label debug
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.printMeshInfo(debug, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    meshRefiner.write(meshRefinement::MESH|meshRefinement::SCALARLEVELS, "");
    if (debug & meshRefinement::OBJINTERSECTIONS)
    {
        meshRefiner.write
        (
            meshRefinement::OBJINTERSECTIONS,
            mesh.time().path()/meshRefiner.timeName()
        );
    }
    Info<< "Wrote mesh in = "
        << mesh.time().cpuTimeIncrement() << " s." << endl;
}



int main(int argc, char *argv[])
{
#   include "addOverwriteOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"

    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    const bool overwrite = args.optionFound("overwrite");

    // Check patches and faceZones are synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    meshRefinement::checkCoupledFaceZones(mesh);


    // Read decomposePar dictionary
    IOdictionary decomposeDict
    (
        IOobject
        (
            "decomposeParDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Read meshing dictionary
    IOdictionary meshDict
    (
       IOobject
       (
            "snappyHexMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    // refinement parameters
    const dictionary& refineDict = meshDict.subDict("castellatedMeshControls");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict.subDict("snapControls");

    // layer addition parameters
    const dictionary& layerDict = meshDict.subDict("addLayersControls");

    // absolute merge distance
    const scalar mergeDist = getMergeDistance
    (
        mesh,
        readScalar(meshDict.lookup("mergeTolerance"))
    );


    // Debug
    // ~~~~~

    const label debug = meshDict.lookupOrDefault<label>("debug", 0);
    if (debug > 0)
    {
        meshRefinement::debug   = debug;
        autoRefineDriver::debug = debug;
        autoSnapDriver::debug   = debug;
        autoLayerDriver::debug  = debug;
    }


    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            //mesh.time().findInstance("triSurface", word::null),// instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict
    );


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement surfaces." << endl;
    refinementSurfaces surfaces
    (
        allGeometry,
        refineDict.subDict("refinementSurfaces")
    );
    Info<< "Read refinement surfaces in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement shells." << endl;
    shellSurfaces shells
    (
        allGeometry,
        refineDict.subDict("refinementRegions")
    );
    Info<< "Read refinement shells in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    Info<< "Setting refinement level of surface to be consistent"
        << " with shells." << endl;
    surfaces.setMinLevelFields(shells);
    Info<< "Checked shell refinement in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Read feature meshes
    // ~~~~~~~~~~~~~~~~~~~

    Info<< "Reading features." << endl;
    refinementFeatures features
    (
        mesh,
        refineDict.lookup("features")
    );
    Info<< "Read features in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;



    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Determining initial surface intersections" << nl
        << "-----------------------------------------" << nl
        << endl;

    // Main refinement engine
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,          // tolerance used in sorting coordinates
        overwrite,          // overwrite mesh files?
        surfaces,           // for surface intersection refinement
        features,           // for feature edges/point based refinement
        shells              // for volume (inside/outside) refinement
    );
    Info<< "Calculated surface intersections in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;

    // Some stats
    meshRefiner.printMeshInfo(debug, "Initial mesh");

    meshRefiner.write
    (
        debug & meshRefinement::OBJINTERSECTIONS,
        mesh.time().path()/meshRefiner.timeName()
    );


    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList globalToPatch;
    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToPatch.setSize(surfaces.nRegions(), -1);

        Info<< "Patch\tType\tRegion" << nl
            << "-----\t----\t------"
            << endl;

        const labelList& surfaceGeometry = surfaces.surfaces();
        const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();

        forAll(surfaceGeometry, surfI)
        {
            label geomI = surfaceGeometry[surfI];

            const wordList& regNames = allGeometry.regionNames()[geomI];

            Info<< surfaces.names()[surfI] << ':' << nl << nl;

            forAll(regNames, i)
            {
                label globalRegionI = surfaces.globalRegion(surfI, i);

                label patchI;

                if (surfacePatchInfo.set(globalRegionI))
                {
                    patchI = meshRefiner.addMeshedPatch
                    (
                        regNames[i],
                        surfacePatchInfo[globalRegionI]
                    );
                }
                else
                {
                    dictionary patchInfo;
                    patchInfo.set("type", wallPolyPatch::typeName);

                    patchI = meshRefiner.addMeshedPatch
                    (
                        regNames[i],
                        patchInfo
                    );
                }

                Info<< patchI << '\t' << mesh.boundaryMesh()[patchI].type()
                    << '\t' << regNames[i] << nl;

                globalToPatch[globalRegionI] = patchI;
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }


    // Parallel
    // ~~~~~~~~

    // Decomposition
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            decomposeDict
        )
    );
    decompositionMethod& decomposer = decomposerPtr();

    if (Pstream::parRun() && !decomposer.parallelAware())
    {
        FatalErrorIn(args.executable())
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware." << endl
            << "Please select one that is (hierarchical, ptscotch)"
            << exit(FatalError);
    }

    // Mesh distribution engine (uses tolerance to reconstruct meshes)
    fvMeshDistribute distributor(mesh, mergeDist);





    // Now do the real work -refinement -snapping -layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const Switch wantRefine(meshDict.lookup("castellatedMesh"));
    const Switch wantSnap(meshDict.lookup("snap"));
    const Switch wantLayers(meshDict.lookup("addLayers"));

    if (wantRefine)
    {
        cpuTime timer;

        autoRefineDriver refineDriver
        (
            meshRefiner,
            decomposer,
            distributor,
            globalToPatch
        );

        // Refinement parameters
        refinementParameters refineParams(refineDict);

        if (!overwrite)
        {
            const_cast<Time&>(mesh.time())++;
        }

        refineDriver.doRefine(refineDict, refineParams, wantSnap, motionDict);

        writeMesh
        (
            "Refined mesh",
            meshRefiner,
            debug
        );

        Info<< "Mesh refined in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantSnap)
    {
        cpuTime timer;

        autoSnapDriver snapDriver
        (
            meshRefiner,
            globalToPatch
        );

        // Snap parameters
        snapParameters snapParams(snapDict);
        // Temporary hack to get access to resolveFeatureAngle
        scalar curvature;
        {
            refinementParameters refineParams(refineDict);
            curvature = refineParams.curvature();
        }

        if (!overwrite)
        {
            const_cast<Time&>(mesh.time())++;
        }

        snapDriver.doSnap(snapDict, motionDict, curvature, snapParams);

        writeMesh
        (
            "Snapped mesh",
            meshRefiner,
            debug
        );

        Info<< "Mesh snapped in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantLayers)
    {
        cpuTime timer;

        autoLayerDriver layerDriver(meshRefiner, globalToPatch);

        // Layer addition parameters
        layerParameters layerParams(layerDict, mesh.boundaryMesh());

        //!!! Temporary hack to get access to maxLocalCells
        bool preBalance;
        {
            refinementParameters refineParams(refineDict);

            preBalance = returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            );
        }


        if (!overwrite)
        {
            const_cast<Time&>(mesh.time())++;
        }

        layerDriver.doLayers
        (
            layerDict,
            motionDict,
            layerParams,
            preBalance,
            decomposer,
            distributor
        );

        writeMesh
        (
            "Layer mesh",
            meshRefiner,
            debug
        );

        Info<< "Layers added in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }


    Info<< "Finished meshing in = "
        << runTime.elapsedCpuTime() << " s." << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
