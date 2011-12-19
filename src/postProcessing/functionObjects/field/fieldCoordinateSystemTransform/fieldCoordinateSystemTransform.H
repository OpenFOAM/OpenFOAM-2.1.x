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
    Foam::fieldCoordinateSystemTransform

Description
    Transforms fields from global cartesian co-ordinates to local co-ordinate
    system

SourceFiles
    fieldCoordinateSystemTransform.C
    IOfieldCoordinateSystemTransform.H

\*---------------------------------------------------------------------------*/

#ifndef fieldCoordinateSystemTransform_H
#define fieldCoordinateSystemTransform_H

#include "OFstream.H"
#include "pointFieldFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
               Class fieldCoordinateSystemTransform Declaration
\*---------------------------------------------------------------------------*/

class fieldCoordinateSystemTransform
{
protected:

    // Protected data

        //- Name
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Fields to transform
        wordList fieldSet_;

        //- Co-ordinate system to transform to
        coordinateSystem coordSys_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        fieldCoordinateSystemTransform(const fieldCoordinateSystemTransform&);

        //- Disallow default bitwise assignment
        void operator=(const fieldCoordinateSystemTransform&);

        template<class Type>
        void transform(const word& fieldName) const;

        template<class Type>
        void transformField(const Type& field) const;


public:

    //- Runtime type information
    TypeName("fieldCoordinateSystemTransform");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fieldCoordinateSystemTransform
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~fieldCoordinateSystemTransform();


    // Member Functions

        //- Return name of the fieldCoordinateSystemTransform object
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the input data
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
#   include "fieldCoordinateSystemTransformTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
