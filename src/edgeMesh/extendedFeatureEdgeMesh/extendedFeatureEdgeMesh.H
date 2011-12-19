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
    Foam::extendedFeatureEdgeMesh

Description

    Description of feature edges and points.

    Feature points are a sorted subset at the start of the overall points list:
        0 .. concaveStart_-1                : convex points (w.r.t normals)
        concaveStart_-1 .. mixedStart_-1    : concave points
        mixedStart_ .. nonFeatureStart_     : mixed internal/external points
        nonFeatureStart_ .. size-1          : non-feature points

    Feature edges are the edgeList of the edgeMesh and are sorted:
        0 .. internalStart_-1           : external edges (convex w.r.t normals)
        internalStart_ .. flatStart_-1  : internal edges (concave)
        flatStart_ .. openStart_-1      : flat edges (neither concave or convex)
                                          can arise from region interfaces on
                                          flat surfaces
        openStart_ .. multipleStart_-1  : open edges (e.g. from baffle surfaces)
        multipleStart_ .. size-1        : multiply connected edges

    The edge direction and feature edge and feature point adjacent normals
    are stored.

SourceFiles
    extendedFeatureEdgeMeshI.H
    extendedFeatureEdgeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFeatureEdgeMesh_H
#define extendedFeatureEdgeMesh_H

#include "edgeMesh.H"
#include "surfaceFeatures.H"
#include "objectRegistry.H"
#include "IOdictionary.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class extendedFeatureEdgeMesh Declaration
\*---------------------------------------------------------------------------*/

class extendedFeatureEdgeMesh
:
    public regIOobject,
    public edgeMesh
{

public:

    //- Runtime type information
    TypeName("extendedFeatureEdgeMesh");

    enum pointStatus
    {
        CONVEX,     // Fully convex point (w.r.t normals)
        CONCAVE,    // Fully concave point
        MIXED,      // A point surrounded by both convex and concave edges
        NONFEATURE  // Not a feature point
    };

    enum edgeStatus
    {
        EXTERNAL,   // "Convex" edge
        INTERNAL,   // "Concave" edge
        FLAT,       // Neither concave or convex, on a flat surface
        OPEN,       // i.e. only connected to one face
        MULTIPLE,   // Multiply connected (connected to more than two faces)
        NONE        // Not a classified feature edge (consistency with
                    // surfaceFeatures)
    };

private:

    // Static data

        //- Angular closeness tolerance for treating normals as the same
        static scalar cosNormalAngleTol_;

        //- Index of the start of the convex feature points - static as 0
        static label convexStart_;

        //- Index of the start of the external feature edges - static as 0
        static label externalStart_;


    // Private data

        //- Index of the start of the concave feature points
        label concaveStart_;

        //- Index of the start of the mixed type feature points
        label mixedStart_;

        //- Index of the start of the non-feature points
        label nonFeatureStart_;

        //- Index of the start of the internal feature edges
        label internalStart_;

        //- Index of the start of the flat feature edges
        label flatStart_;

        //- Index of the start of the open feature edges
        label openStart_;

        //- Index of the start of the multiply-connected feature edges
        label multipleStart_;

        //- Normals of the features, to be referred to by index by both feature
        //  points and edges, unsorted
        vectorField normals_;

        //- Flat and open edges require the direction of the edge
        vectorField edgeDirections_;

        //- Indices of the normals that are adjacent to the feature edges
        labelListList edgeNormals_;

        //- Indices of the normals that are adjacent to the feature points
        labelListList featurePointNormals_;

        //- Feature edges which are on the boundary between regions
        labelList regionEdges_;

        //- Search tree for all edges
        mutable autoPtr<indexedOctree<treeDataEdge> > edgeTree_;

        //- Individual search trees for each type of edge
        mutable PtrList<indexedOctree<treeDataEdge> > edgeTreesByType_;


    // Private Member Functions

        //- Classify the type of feature point.  Requires valid stored member
        //  data for edges and normals.
        pointStatus classifyFeaturePoint(label ptI) const;

        //- Classify the type of feature edge.  Requires face centre 0 to face
        //  centre 1 vector to distinguish internal from external
        edgeStatus classifyEdge
        (
            const List<vector>& norms,
            const labelList& edNorms,
            const vector& fC0tofC1
        ) const;


public:

    // Static data

        //- Number of possible point types (i.e. number of slices)
        static label nPointTypes;

        //- Number of possible feature edge types (i.e. number of slices)
        static label nEdgeTypes;

    // Constructors

        //- Construct (read) given an IOobject
        extendedFeatureEdgeMesh(const IOobject&);

        //- Construct as copy
        explicit extendedFeatureEdgeMesh
        (
            const IOobject&,
            const extendedFeatureEdgeMesh&
        );

        //- Construct by transferring components (points, edges)
        extendedFeatureEdgeMesh
        (
            const IOobject&,
            const Xfer<pointField>&,
            const Xfer<edgeList>&
        );

        //- Construct (read) given surfaceFeatures, an objectRegistry and a
        //  fileName to write to.  Extracts, classifies and reorders the data
        //  from surfaceFeatures.
        extendedFeatureEdgeMesh
        (
            const surfaceFeatures& sFeat,
            const objectRegistry& obr,
            const fileName& sFeatFileName
        );

        //- Construct from all components
        extendedFeatureEdgeMesh
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
        );


    //- Destructor
    ~extendedFeatureEdgeMesh();


    // Member Functions

        // Find

            //- Find nearest surface edge for the sample point.
            void nearestFeatureEdge
            (
                const point& sample,
                scalar searchDistSqr,
                pointIndexHit& info
            ) const;

            //- Find nearest surface edge for each sample point.
            void nearestFeatureEdge
            (
                const pointField& samples,
                const scalarField& searchDistSqr,
                List<pointIndexHit>& info
            ) const;

            //- Find the nearest point on each type of feature edge
            void nearestFeatureEdgeByType
            (
                const point& sample,
                const scalarField& searchDistSqr,
                List<pointIndexHit>& info
            ) const;

        // Access

            //- Return the index of the start of the convex feature points
            inline label convexStart() const;

            //- Return the index of the start of the concave feature points
            inline label concaveStart() const;

            //- Return the index of the start of the mixed type feature points
            inline label mixedStart() const;

            //- Return the index of the start of the non-feature points
            inline label nonFeatureStart() const;

            //- Return the index of the start of the external feature edges
            inline label externalStart() const;

            //- Return the index of the start of the internal feature edges
            inline label internalStart() const;

            //- Return the index of the start of the flat feature edges
            inline label flatStart() const;

            //- Return the index of the start of the open feature edges
            inline label openStart() const;

            //- Return the index of the start of the multiply-connected feature
            //  edges
            inline label multipleStart() const;

            //- Return whether or not the point index is a feature point
            inline bool featurePoint(label ptI) const;

            //- Return the normals of the surfaces adjacent to the feature edges
            //  and points
            inline const vectorField& normals() const;

            //- Return the edgeDirection vectors
            inline const vectorField& edgeDirections() const;

            //- Return the direction of edgeI, pointing away from ptI
            inline vector edgeDirection(label edgeI, label ptI) const;

            //- Return the indices of the normals that are adjacent to the
            //  feature edges
            inline const labelListList& edgeNormals() const;

            //- Return the normal vectors for a given set of normal indices
            inline vectorField edgeNormals(const labelList& edgeNormIs) const;

            //- Return the normal vectors for a given edge
            inline vectorField edgeNormals(label edgeI) const;

            //- Return the indices of the normals that are adjacent to the
            //  feature points
            inline const labelListList& featurePointNormals() const;

            //- Return the normal vectors for a given feature point
            inline vectorField featurePointNormals(label ptI) const;

            //- Return the feature edges which are on the boundary between
            //  regions
            inline const labelList& regionEdges() const;

            //- Return the pointStatus of a specified point
            inline pointStatus getPointStatus(label ptI) const;

            //- Return the edgeStatus of a specified edge
            inline edgeStatus getEdgeStatus(label edgeI) const;

            //- Demand driven construction of octree for boundary edges
            const indexedOctree<treeDataEdge>& edgeTree() const;

            //- Demand driven construction of octree for boundary edges by type
            const PtrList<indexedOctree<treeDataEdge> >&
            edgeTreesByType() const;


        // Write

            //- Write all components of the extendedFeatureEdgeMesh as obj files
            void writeObj(const fileName& prefix) const;

            //- Give precedence to the regIOobject write
            using regIOobject::write;

            //- WriteData function required for regIOobject write operation
            virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "extendedFeatureEdgeMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
