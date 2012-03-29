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

Description
    Makes internal faces into boundary faces. Does not duplicate points, unlike
    mergeOrSplitBaffles.

    Note: if any coupled patch face is selected for baffling the opposite
    member has to be selected for baffling as well. Note that this
    is the same as repatching. This was added only for convenience so
    you don't have to filter coupled boundary out of your set.

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "Time.H"
#include "faceSet.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ZoneIDs.H"
#include "fvMeshMapper.H"
#include "SetPatchFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void modifyOrAddFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label faceI,
    const label own,
    const bool flipFaceFlux,
    const label newPatchI,
    const label zoneID,
    const bool zoneFlip,

    PackedBoolList& modifiedFace
)
{
    if (!modifiedFace[faceI])
    {
        // First usage of face. Modify.
        meshMod.setAction
        (
            polyModifyFace
            (
                f,                          // modified face
                faceI,                      // label of face
                own,                        // owner
                -1,                         // neighbour
                flipFaceFlux,               // face flip
                newPatchI,                  // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
        modifiedFace[faceI] = 1;
    }
    else
    {
        // Second or more usage of face. Add.
        meshMod.setAction
        (
            polyAddFace
            (
                f,                          // modified face
                own,                        // owner
                -1,                         // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face
                flipFaceFlux,               // face flip
                newPatchI,                  // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
}


label findPatchID(const polyMesh& mesh, const word& name)
{
    const label patchI = mesh.boundaryMesh().findPatchID(name);

    if (patchI == -1)
    {
        FatalErrorIn("findPatchID(const polyMesh&, const word&)")
            << "Cannot find patch " << name << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }
    return patchI;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Makes internal faces into boundary faces.\n"
        "Does not duplicate points, unlike mergeOrSplitBaffles."
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"

    argList::validArgs.append("faceZone");
    argList::validArgs.append("(masterPatch slavePatch)");
    argList::addOption
    (
        "additionalPatches",
        "((master2 slave2) .. (masterN slaveN))"
    );
    argList::addBoolOption
    (
        "internalFacesOnly",
        "do not convert boundary faces"
    );
    argList::addBoolOption
    (
        "updateFields",
        "update fields to include new patches:"
        " NOTE: updated field values may need to be edited"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const faceZoneMesh& faceZones = mesh.faceZones();

    // Faces to baffle
    faceZoneID zoneID(args.additionalArgs()[0], faceZones);

    Info<< "Converting faces on zone " << zoneID.name()
        << " into baffles." << nl << endl;

    if (zoneID.index() == -1)
    {
        FatalErrorIn(args.executable()) << "Cannot find faceZone "
            << zoneID.name() << endl
            << "Valid zones are " << faceZones.names()
            << exit(FatalError);
    }

    const faceZone& fZone = faceZones[zoneID.index()];

    Info<< "Found " << returnReduce(fZone.size(), sumOp<label>())
        << " faces on zone " << zoneID.name() << nl << endl;

    // Make sure patches and zoneFaces are synchronised across couples
    patches.checkParallelSync(true);
    fZone.checkParallelSync(true);

    // Patches to put baffles into
    DynamicList<label> newMasterPatches(1);
    DynamicList<label> newSlavePatches(1);

    const Pair<word> patchNames(IStringStream(args.additionalArgs()[1])());
    newMasterPatches.append(findPatchID(mesh, patchNames[0]));
    newSlavePatches.append(findPatchID(mesh, patchNames[1]));
    Info<< "Using master patch " << patchNames[0]
        << " at index " << newMasterPatches[0] << endl;
    Info<< "Using slave patch " << patchNames[1]
        << " at index " << newSlavePatches[0] << endl;


    // Additional patches
    if (args.optionFound("additionalPatches"))
    {
        const List<Pair<word> > patchNames
        (
            args.optionLookup("additionalPatches")()
        );

        newMasterPatches.reserve(patchNames.size() + 1);
        newSlavePatches.reserve(patchNames.size() + 1);
        forAll(patchNames, i)
        {
            newMasterPatches.append(findPatchID(mesh, patchNames[i][0]));
            newSlavePatches.append(findPatchID(mesh, patchNames[i][1]));
            Info<< "Using additional patches " << patchNames[i]
                << " at indices " << newMasterPatches.last()
                << " and " << newSlavePatches.last()
                << endl;
        }
    }


    const bool overwrite = args.optionFound("overwrite");
    const bool internalFacesOnly = args.optionFound("internalFacesOnly");

    if (internalFacesOnly)
    {
        Info<< "Not converting faces on non-coupled patches." << nl << endl;
    }


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.
    Info<< "Reading geometric fields" << nl << endl;
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);


    // Mesh change container
    polyTopoChange meshMod(mesh);


    // Do the actual changes. Note:
    // - loop in incrementing face order (not necessary if faceZone ordered).
    //   Preserves any existing ordering on patch faces.
    // - two passes, do non-flip faces first and flip faces second. This
    //   guarantees that when e.g. creating a cyclic all faces from one
    //   side come first and faces from the other side next.

    // Whether first use of face (modify) or consecutive (add)
    PackedBoolList modifiedFace(mesh.nFaces());
    label nModified = 0;

    forAll(newMasterPatches, i)
    {
        // Pass 1. Do selected side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            label zoneFaceI = fZone.whichFace(faceI);

            if (zoneFaceI != -1)
            {
                if (!fZone.flipMap()[zoneFaceI])
                {
                    // Use owner side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[faceI],    // modified face
                        faceI,                  // label of face
                        mesh.faceOwner()[faceI],// owner
                        false,                  // face flip
                        newMasterPatches[i],    // patch for face
                        zoneID.index(),         // zone for face
                        false,                  // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
                else
                {
                    // Use neighbour side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[faceI].reverseFace(),  // modified face
                        faceI,                      // label of face
                        mesh.faceNeighbour()[faceI],// owner
                        true,                       // face flip
                        newMasterPatches[i],        // patch for face
                        zoneID.index(),             // zone for face
                        true,                       // face flip in zone
                        modifiedFace                // modify or add status
                    );
                }

                nModified++;
            }
        }


        // Pass 2. Do other side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            label zoneFaceI = fZone.whichFace(faceI);

            if (zoneFaceI != -1)
            {
                if (!fZone.flipMap()[zoneFaceI])
                {
                    // Use neighbour side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[faceI].reverseFace(),  // modified face
                        faceI,                              // label of face
                        mesh.faceNeighbour()[faceI],        // owner
                        true,                               // face flip
                        newSlavePatches[i],                 // patch for face
                        zoneID.index(),                     // zone for face
                        true,                               // face flip in zone
                        modifiedFace                        // modify or add
                    );
                }
                else
                {
                    // Use owner side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[faceI],    // modified face
                        faceI,                  // label of face
                        mesh.faceOwner()[faceI],// owner
                        false,                  // face flip
                        newSlavePatches[i],     // patch for face
                        zoneID.index(),         // zone for face
                        false,                  // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
            }
        }


        // Modify any boundary faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        // Normal boundary:
        // - move to new patch. Might already be back-to-back baffle
        // you want to add cyclic to. Do warn though.
        //
        // Processor boundary:
        // - do not move to cyclic
        // - add normal patches though.

        // For warning once per patch.
        labelHashSet patchWarned;

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            label newPatchI = newMasterPatches[i];

            if (pp.coupled() && patches[newPatchI].coupled())
            {
                // Do not allow coupled faces to be moved to different coupled
                // patches.
            }
            else if (pp.coupled() || !internalFacesOnly)
            {
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;

                    label zoneFaceI = fZone.whichFace(faceI);

                    if (zoneFaceI != -1)
                    {
                        if (patchWarned.insert(patchI))
                        {
                            WarningIn(args.executable())
                                << "Found boundary face (in patch " << pp.name()
                                << ") in faceZone " << fZone.name()
                                << " to convert to baffle patch "
                                << patches[newPatchI].name()
                                << endl
                                << "    Run with -internalFacesOnly option"
                                << " if you don't wish to convert"
                                << " boundary faces." << endl;
                        }

                        modifyOrAddFace
                        (
                            meshMod,
                            mesh.faces()[faceI],        // modified face
                            faceI,                      // label of face
                            mesh.faceOwner()[faceI],    // owner
                            false,                      // face flip
                            newPatchI,                  // patch for face
                            zoneID.index(),             // zone for face
                            fZone.flipMap()[zoneFaceI], // face flip in zone
                            modifiedFace                // modify or add status
                        );
                        nModified++;
                    }
                }
            }
        }
    }


    Info<< "Converted " << returnReduce(nModified, sumOp<label>())
        << " faces into boundary faces on patches " << patchNames << nl << endl;

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Correct boundary faces mapped-out-of-nothing.
    {
        fvMeshMapper mapper(mesh, map);
        bool hasWarned = false;
        forAll(newMasterPatches, i)
        {
            label patchI = newMasterPatches[i];
            const fvPatchMapper& pm = mapper.boundaryMap()[patchI];
            if (pm.sizeBeforeMapping() == 0)
            {
                if (!hasWarned)
                {
                    hasWarned = true;
                    WarningIn(args.executable())
                        << "Setting field on boundary faces to zero." << endl
                        << "You might have to edit these fields." << endl;
                }

                SetPatchFields(vsFlds, patchI, pTraits<scalar>::zero);
                SetPatchFields(vvFlds, patchI, pTraits<vector>::zero);
                SetPatchFields(vstFlds, patchI, pTraits<sphericalTensor>::zero);
                SetPatchFields(vsymtFlds, patchI, pTraits<symmTensor>::zero);
                SetPatchFields(vtFlds, patchI, pTraits<tensor>::zero);

                SetPatchFields(ssFlds, patchI, pTraits<scalar>::zero);
                SetPatchFields(svFlds, patchI, pTraits<vector>::zero);
                SetPatchFields(sstFlds, patchI, pTraits<sphericalTensor>::zero);
                SetPatchFields(ssymtFlds, patchI, pTraits<symmTensor>::zero);
                SetPatchFields(stFlds, patchI, pTraits<tensor>::zero);
            }
        }
        forAll(newSlavePatches, i)
        {
            label patchI = newSlavePatches[i];
            const fvPatchMapper& pm = mapper.boundaryMap()[patchI];
            if (pm.sizeBeforeMapping() == 0)
            {
                SetPatchFields(vsFlds, patchI, pTraits<scalar>::zero);
                SetPatchFields(vvFlds, patchI, pTraits<vector>::zero);
                SetPatchFields(vstFlds, patchI, pTraits<sphericalTensor>::zero);
                SetPatchFields(vsymtFlds, patchI, pTraits<symmTensor>::zero);
                SetPatchFields(vtFlds, patchI, pTraits<tensor>::zero);

                SetPatchFields(ssFlds, patchI, pTraits<scalar>::zero);
                SetPatchFields(svFlds, patchI, pTraits<vector>::zero);
                SetPatchFields(sstFlds, patchI, pTraits<sphericalTensor>::zero);
                SetPatchFields(ssymtFlds, patchI, pTraits<symmTensor>::zero);
                SetPatchFields(stFlds, patchI, pTraits<tensor>::zero);
            }
        }
    }

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    Info<< "Writing mesh to " << runTime.timeName() << endl;

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
