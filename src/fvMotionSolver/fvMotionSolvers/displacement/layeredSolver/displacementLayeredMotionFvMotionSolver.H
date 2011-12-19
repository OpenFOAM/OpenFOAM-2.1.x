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

Class
    Foam::displacementLayeredMotionFvMotionSolver

Description
    Mesh motion solver for an (multi-block) extruded fvMesh. Gets given the
    structure of the mesh blocks and boundary conditions on these blocks.

    Note: should not be an fvMotionSolver but just a motionSolver. Only here
    so we can reuse displacementFvMotionSolver functionality (e.g. surface
    following boundary conditions)


    The displacementLayeredMotionCoeffs subdict of dynamicMeshDict specifies
    per region (=cellZone) the boundary conditions on two opposing patches
    (=faceZones). It then interpolates the boundary values using topological
    walking so requires the cellZone to be a layered mesh.

    The boundary conditions on faceZones are currently:

    follow: the faceZone follows the overall mesh displacement.
            Use this for faceZones on boundary faces (so it uses the
            proper boundary conditions on the pointDisplacement).

    fixedValue: fixed value.

    timeVaryingUniformFixedValue: table-driven fixed value.

    slip: the second of a pair of faceZones follows the tangential movement
          specified by the first faceZone. (= removes the normal component).

SourceFiles
    displacementLayeredMotionFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementLayeredMotionFvMotionSolver_H
#define displacementLayeredMotionFvMotionSolver_H

#include "displacementFvMotionSolver.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations

/*---------------------------------------------------------------------------*\
             Class displacementLayeredMotionFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementLayeredMotionFvMotionSolver
:
    public displacementFvMotionSolver
{
    // Private data

        //- Point motion field
        mutable pointVectorField pointDisplacement_;


    // Private Member Functions

        void calcZoneMask
        (
            const label cellZoneI,
            PackedBoolList& isZonePoint,
            PackedBoolList& isZoneEdge
        ) const;

        void walkStructured
        (
            const label cellZoneI,
            const PackedBoolList& isZonePoint,
            const PackedBoolList& isZoneEdge,
            const labelList& seedPoints,
            const vectorField& seedData,
            scalarField& distance,
            vectorField& data
        ) const;

        tmp<vectorField> faceZoneEvaluate
        (
            const faceZone& fz,
            const labelList& meshPoints,
            const dictionary& dict,
            const PtrList<pointVectorField>& patchDisp,
            const label patchI
        ) const;

        void cellZoneSolve
        (
            const label cellZoneI,
            const dictionary& zoneDict
        );


        //- Disallow default bitwise copy construct
        displacementLayeredMotionFvMotionSolver
        (
            const displacementLayeredMotionFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const displacementLayeredMotionFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("displacementLayeredMotion");


    // Constructors

        //- Construct from polyMesh and data stream
        displacementLayeredMotionFvMotionSolver
        (
            const polyMesh&,
            Istream& msDataUnused
        );


    //- Destructor
    ~displacementLayeredMotionFvMotionSolver();


    // Member Functions

        //- Return reference to the point motion displacement field
        pointVectorField& pointDisplacement()
        {
            return pointDisplacement_;
        }

        //- Return const reference to the point motion displacement field
        const pointVectorField& pointDisplacement() const
        {
            return pointDisplacement_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
