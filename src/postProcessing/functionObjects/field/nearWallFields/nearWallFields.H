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
    Foam::nearWallFields

Description
    Samples near-patch volFields

    Holds fields
    - every timestep the field get updated with new values
    - at write it writes the fields
    so this functionObject can either be used to calculate a new field
    as a postprocessing step or (since the fields are registered)
    use these in another functionObject (e.g. faceSource).

    surfaceValues
    {
        type        nearWallFields;
        ..
        enabled         true;
        outputControl   outputTime;
        ..
        // Name of volField and corresponding surfaceField
        fields      ((p pNear)(U UNear));
        // Name of patch to sample
        patches     (movingWall);
        // Distance away from the wall
        distance    0.13;   // distance away from wall
    }


SourceFiles
    nearWallFields.C
    IOnearWallFields.H

\*---------------------------------------------------------------------------*/

#ifndef nearWallFields_H
#define nearWallFields_H

#include "OFstream.H"
#include "volFields.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class nearWallFields Declaration
\*---------------------------------------------------------------------------*/

class nearWallFields
{
protected:

    // Protected data

        //- Name of this set of nearWallFields object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        // Read from dictionary

            //- Fields to process
            List<Tuple2<word, word> > fieldSet_;

            //- Patches to sample
            labelHashSet patchSet_;

            //- Distance away from wall
            scalar distance_;

            //- From original field to sampled result
            HashTable<word> fieldMap_;

            //- From resulting back to original field
            HashTable<word> reverseFieldMap_;

        //- Locally constructed fields
        PtrList<volScalarField> vsf_;
        PtrList<volVectorField> vvf_;
        PtrList<volSphericalTensorField> vSpheretf_;
        PtrList<volSymmTensorField> vSymmtf_;
        PtrList<volTensorField> vtf_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        nearWallFields(const nearWallFields&);

        //- Disallow default bitwise assignment
        void operator=(const nearWallFields&);

        template<class Type>
        void createFields
        (
            PtrList<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;

        template<class Type>
        void sampleFields
        (
            PtrList<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;


public:

    //- Runtime type information
    TypeName("nearWallFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        nearWallFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~nearWallFields();


    // Member Functions

        //- Return name of the nearWallFields object
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "nearWallFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
