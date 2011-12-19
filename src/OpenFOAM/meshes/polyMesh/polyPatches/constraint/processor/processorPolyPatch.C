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

#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "matchPoints.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "Time.H"
#include "transformList.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo
)
:
    coupledPolyPatch(name, size, start, index, bm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::~processorPolyPatch()
{
    neighbPointsPtr_.clear();
    neighbEdgesPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << faceCentres()
            << faceAreas()
            << faceCellCentres();
    }
}


void Foam::processorPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        {
            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> neighbFaceCentres_
                >> neighbFaceAreas_
                >> neighbFaceCellCentres_;
        }


        // My normals
        vectorField faceNormals(size());

        // Neighbour normals
        vectorField nbrFaceNormals(neighbFaceAreas_.size());

        // Face match tolerances
        scalarField tols =
            calcFaceTol(*this, points(), faceCentres());

        // Calculate normals from areas and check
        forAll(faceNormals, facei)
        {
            scalar magSf = mag(faceAreas()[facei]);
            scalar nbrMagSf = mag(neighbFaceAreas_[facei]);
            scalar avSf = (magSf + nbrMagSf)/2.0;

            if (magSf < ROOTVSMALL && nbrMagSf < ROOTVSMALL)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check. (note use of sqrt(VSMALL) since that is how mag
                // scales)
                faceNormals[facei] = point(1, 0, 0);
                nbrFaceNormals[facei] = faceNormals[facei];
            }
            else if (mag(magSf - nbrMagSf) > matchTolerance()*sqr(tols[facei]))
            {
                fileName nm
                (
                    boundaryMesh().mesh().time().path()
                   /name()+"_faces.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry : Writing my "
                    << size()
                    << " faces to OBJ file " << nm << endl;

                writeOBJ(nm, *this, points());

                OFstream ccStr
                (
                    boundaryMesh().mesh().time().path()
                    /name() + "_faceCentresConnections.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry :"
                    << " Dumping cell centre lines between"
                    << " corresponding face centres to OBJ file" << ccStr.name()
                    << endl;

                label vertI = 0;

                forAll(faceCentres(), faceI)
                {
                    const point& c0 = neighbFaceCentres_[faceI];
                    const point& c1 = faceCentres()[faceI];

                    writeOBJ(ccStr, c0, c1, vertI);
                }

                FatalErrorIn
                (
                    "processorPolyPatch::calcGeometry()"
                )   << "face " << facei << " area does not match neighbour by "
                    << 100*mag(magSf - nbrMagSf)/avSf
                    << "% -- possible face ordering problem." << endl
                    << "patch:" << name()
                    << " my area:" << magSf
                    << " neighbour area:" << nbrMagSf
                    << " matching tolerance:"
                    << matchTolerance()*sqr(tols[facei])
                    << endl
                    << "Mesh face:" << start()+facei
                    << " vertices:"
                    << UIndirectList<point>(points(), operator[](facei))()
                    << endl
                    << "If you are certain your matching is correct"
                    << " you can increase the 'matchTolerance' setting"
                    << " in the patch dictionary in the boundary file."
                    << endl
                    << "Rerun with processor debug flag set for"
                    << " more information." << exit(FatalError);
            }
            else
            {
                faceNormals[facei] = faceAreas()[facei]/magSf;
                nbrFaceNormals[facei] = neighbFaceAreas_[facei]/nbrMagSf;
            }
        }

        calcTransformTensors
        (
            faceCentres(),
            neighbFaceCentres_,
            faceNormals,
            nbrFaceNormals,
            matchTolerance()*tols,
            matchTolerance()
        );
    }
}


void Foam::processorPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    processorPolyPatch::initGeometry(pBufs);
}


void Foam::processorPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField&
)
{
    processorPolyPatch::calcGeometry(pBufs);
}


void Foam::processorPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);

    if (Pstream::parRun())
    {
        // Express all points as patch face and index in face.
        labelList pointFace(nPoints());
        labelList pointIndex(nPoints());

        for (label patchPointI = 0; patchPointI < nPoints(); patchPointI++)
        {
            label faceI = pointFaces()[patchPointI][0];

            pointFace[patchPointI] = faceI;

            const face& f = localFaces()[faceI];

            pointIndex[patchPointI] = findIndex(f, patchPointI);
        }

        // Express all edges as patch face and index in face.
        labelList edgeFace(nEdges());
        labelList edgeIndex(nEdges());

        for (label patchEdgeI = 0; patchEdgeI < nEdges(); patchEdgeI++)
        {
            label faceI = edgeFaces()[patchEdgeI][0];

            edgeFace[patchEdgeI] = faceI;

            const labelList& fEdges = faceEdges()[faceI];

            edgeIndex[patchEdgeI] = findIndex(fEdges, patchEdgeI);
        }

        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << pointFace
            << pointIndex
            << edgeFace
            << edgeIndex;
    }
}


void Foam::processorPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    // For completeness
    polyPatch::updateMesh(pBufs);

    neighbPointsPtr_.clear();
    neighbEdgesPtr_.clear();

    if (Pstream::parRun())
    {
        labelList nbrPointFace;
        labelList nbrPointIndex;
        labelList nbrEdgeFace;
        labelList nbrEdgeIndex;

        {
            // Note cannot predict exact size since opposite nPoints might
            // be different from one over here.
            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> nbrPointFace
                >> nbrPointIndex
                >> nbrEdgeFace
                >> nbrEdgeIndex;
        }

        // Convert neighbour faces and indices into face back into
        // my edges and points.

        // Convert points.
        // ~~~~~~~~~~~~~~~

        neighbPointsPtr_.reset(new labelList(nPoints(), -1));
        labelList& neighbPoints = neighbPointsPtr_();

        forAll(nbrPointFace, nbrPointI)
        {
            // Find face and index in face on this side.
            const face& f = localFaces()[nbrPointFace[nbrPointI]];
            label index = (f.size() - nbrPointIndex[nbrPointI]) % f.size();
            label patchPointI = f[index];

            if (neighbPoints[patchPointI] == -1)
            {
                // First reference of point
                neighbPoints[patchPointI] = nbrPointI;
            }
            else if (neighbPoints[patchPointI] >= 0)
            {
                // Point already visited. Mark as duplicate.
                neighbPoints[patchPointI] = -2;
            }
        }

        // Reset all duplicate entries to -1.
        forAll(neighbPoints, patchPointI)
        {
            if (neighbPoints[patchPointI] == -2)
            {
                neighbPoints[patchPointI] = -1;
            }
        }

        // Convert edges.
        // ~~~~~~~~~~~~~~

        neighbEdgesPtr_.reset(new labelList(nEdges(), -1));
        labelList& neighbEdges = neighbEdgesPtr_();

        forAll(nbrEdgeFace, nbrEdgeI)
        {
            // Find face and index in face on this side.
            const labelList& f = faceEdges()[nbrEdgeFace[nbrEdgeI]];
            label index = (f.size() - nbrEdgeIndex[nbrEdgeI] - 1) % f.size();
            label patchEdgeI = f[index];

            if (neighbEdges[patchEdgeI] == -1)
            {
                // First reference of edge
                neighbEdges[patchEdgeI] = nbrEdgeI;
            }
            else if (neighbEdges[patchEdgeI] >= 0)
            {
                // Edge already visited. Mark as duplicate.
                neighbEdges[patchEdgeI] = -2;
            }
        }

        // Reset all duplicate entries to -1.
        forAll(neighbEdges, patchEdgeI)
        {
            if (neighbEdges[patchEdgeI] == -2)
            {
                neighbEdges[patchEdgeI] = -1;
            }
        }

        // Remove any addressing used for shared points/edges calculation
        // since mostly not needed.
        primitivePatch::clearOut();
    }
}


const Foam::labelList& Foam::processorPolyPatch::neighbPoints() const
{
    if (!neighbPointsPtr_.valid())
    {
        FatalErrorIn("processorPolyPatch::neighbPoints() const")
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return neighbPointsPtr_();
}


const Foam::labelList& Foam::processorPolyPatch::neighbEdges() const
{
    if (!neighbEdgesPtr_.valid())
    {
        FatalErrorIn("processorPolyPatch::neighbEdges() const")
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return neighbEdgesPtr_();
}


void Foam::processorPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        fileName nm
        (
            boundaryMesh().mesh().time().path()
           /name()+"_faces.obj"
        );
        Pout<< "processorPolyPatch::order : Writing my " << pp.size()
            << " faces to OBJ file " << nm << endl;
        writeOBJ(nm, pp, pp.points());

        // Calculate my face centres
        const pointField& fc = pp.faceCentres();

        OFstream localStr
        (
            boundaryMesh().mesh().time().path()
           /name() + "_localFaceCentres.obj"
        );
        Pout<< "processorPolyPatch::order : "
            << "Dumping " << fc.size()
            << " local faceCentres to " << localStr.name() << endl;

        forAll(fc, faceI)
        {
            writeOBJ(localStr, fc[faceI]);
        }
    }

    if (owner())
    {
        pointField anchors(getAnchorPoints(pp, pp.points()));

        // Now send all info over to the neighbour
        UOPstream toNeighbour(neighbProcNo(), pBufs);
        toNeighbour << pp.faceCentres() << anchors;
    }
}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool Foam::processorPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    // Note: we only get the faces that originate from internal faces.

    if (!Pstream::parRun())
    {
        return false;
    }

    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (owner())
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFaceI)
        {
            faceMap[patchFaceI] = patchFaceI;
        }

        return false;
    }
    else
    {
        vectorField masterCtrs;
        vectorField masterAnchors;

        // Receive data from neighbour
        {
            UIPstream fromNeighbour(neighbProcNo(), pBufs);
            fromNeighbour >> masterCtrs >> masterAnchors;
        }

        // Calculate typical distance from face centre
        scalarField tols
        (
            matchTolerance()*calcFaceTol(pp, pp.points(), pp.faceCentres())
        );

        if (debug || masterCtrs.size() != pp.size())
        {
            {
                OFstream nbrStr
                (
                    boundaryMesh().mesh().time().path()
                   /name() + "_nbrFaceCentres.obj"
                );
                Pout<< "processorPolyPatch::order : "
                    << "Dumping neighbour faceCentres to " << nbrStr.name()
                    << endl;
                forAll(masterCtrs, faceI)
                {
                    writeOBJ(nbrStr, masterCtrs[faceI]);
                }
            }

            if (masterCtrs.size() != pp.size())
            {
                FatalErrorIn
                (
                    "processorPolyPatch::order(const primitivePatch&"
                    ", labelList&, labelList&) const"
                )   << "in patch:" << name() << " : "
                    << "Local size of patch is " << pp.size() << " (faces)."
                    << endl
                    << "Received from neighbour " << masterCtrs.size()
                    << " faceCentres!"
                    << abort(FatalError);
            }
        }

        // Geometric match of face centre vectors
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // 1. Try existing ordering and transformation
        bool matchedAll = matchPoints
        (
            pp.faceCentres(),
            masterCtrs,
            tols,
            true,
            faceMap
        );

        if (!matchedAll || debug)
        {
            // Dump faces
            fileName str
            (
                boundaryMesh().mesh().time().path()
               /name()/name()+"_faces.obj"
            );
            Pout<< "processorPolyPatch::order :"
                << " Writing faces to OBJ file " << str.name() << endl;
            writeOBJ(str, pp, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /name() + "_faceCentresConnections.obj"
            );

            Pout<< "processorPolyPatch::order :"
                << " Dumping newly found match as lines between"
                << " corresponding face centres to OBJ file " << ccStr.name()
                << endl;

            label vertI = 0;

            forAll(pp.faceCentres(), faceI)
            {
                label masterFaceI = faceMap[faceI];

                if (masterFaceI != -1)
                {
                    const point& c0 = masterCtrs[masterFaceI];
                    const point& c1 = pp.faceCentres()[faceI];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }

        if (!matchedAll)
        {
            SeriousErrorIn
            (
                "processorPolyPatch::order(const primitivePatch&"
                ", labelList&, labelList&) const"
            )   << "in patch:" << name() << " : "
                << "Cannot match vectors to faces on both sides of patch"
                << endl
                << "    masterCtrs[0]:" << masterCtrs[0] << endl
                << "    ctrs[0]:" << pp.faceCentres()[0] << endl
                << "    Please check your topology changes or maybe you have"
                << " multiple separated (from cyclics) processor patches"
                << endl
                << "    Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }

        // Set rotation.
        forAll(faceMap, oldFaceI)
        {
            // The face f will be at newFaceI (after morphing) and we want its
            // anchorPoint (= f[0]) to align with the anchorpoint for the
            // corresponding face on the other side.

            label newFaceI = faceMap[oldFaceI];

            const point& wantedAnchor = masterAnchors[newFaceI];

            rotation[newFaceI] = getRotation
            (
                pp.points(),
                pp[oldFaceI],
                wantedAnchor,
                tols[oldFaceI]
            );

            if (rotation[newFaceI] == -1)
            {
                SeriousErrorIn
                (
                    "processorPolyPatch::order(const primitivePatch&"
                    ", labelList&, labelList&) const"
                )   << "in patch " << name()
                    << " : "
                    << "Cannot find point on face " << pp[oldFaceI]
                    << " with vertices "
                    << UIndirectList<point>(pp.points(), pp[oldFaceI])()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of processor patch " << name()
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;

                return false;
            }
        }

        forAll(faceMap, faceI)
        {
            if (faceMap[faceI] != faceI || rotation[faceI] != 0)
            {
                return true;
            }
        }

        return false;
    }
}


void Foam::processorPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    os.writeKeyword("myProcNo") << myProcNo_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbProcNo") << neighbProcNo_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
