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
    Foam::regionModel

Description
    Base class for region models

SourceFiles
    regionModelI.H
    regionModel.C

\*---------------------------------------------------------------------------*/

#ifndef regionModel_H
#define regionModel_H

#include "IOdictionary.H"
#include "Switch.H"
#include "labelList.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
//class fvMesh;
//class Time;

namespace regionModels
{

/*---------------------------------------------------------------------------*\
                       Class regionModel Declaration
\*---------------------------------------------------------------------------*/

class regionModel
:
    public IOdictionary
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionModel(const regionModel&);

        //- Disallow default bitwise assignment
        void operator=(const regionModel&);

        //- Construct region mesh and fields
        void constructMeshObjects();

        //- Construct region mesh and dictionary
        void constructMeshObjects(const dictionary& dict);

        //- Initialise the region
        void initialise();


protected:

    // Protected data

        //- Reference to the primary mesh database
        const fvMesh& primaryMesh_;

        //- Reference to the time database
        const Time& time_;

        //- Active flag
        Switch active_;

        //- Active information output
        Switch infoOutput_;

        //- Model name
        const word modelName_;

        //- Pointer to the region mesh database
        autoPtr<fvMesh> regionMeshPtr_;

        //- Model coefficients dictionary
        dictionary coeffs_;


        // Addressing

            //- List of patch IDs on the primary region coupled to this region
            labelList primaryPatchIDs_;

            //- List of patch IDs internally coupled with the primary region
            labelList intCoupledPatchIDs_;


        //- Region name
        word regionName_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("regionModel");


    // Constructors

        //- Construct null
        regionModel(const fvMesh& mesh);

        //- Construct from mesh, region type and name
        regionModel
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            bool readFields = true
        );

        //- Construct from mesh and name and dict
        regionModel
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            const dictionary& dict,
            bool readFields = true
        );



    //- Destructor
    virtual ~regionModel();


    // Member Functions

        // Access

            //- Return the reference to the primary mesh database
            inline const fvMesh& primaryMesh() const;

            //- Return the reference to the time database
            inline const Time& time() const;

            //- Return the active flag
            inline const Switch& active() const;

            //- Return the information flag
            inline const Switch& infoOutput() const;

            //- Return the model name
            inline const word& modelName() const;

            //- Return the region mesh database
            inline const fvMesh& regionMesh() const;

            //- Return the region mesh database for manipulation
            inline fvMesh& regionMesh();

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;

            //- Return the solution dictionary
            inline const dictionary& solution() const;


            // Addressing

                //- Return true if patchI on the primary region is a coupled
                //  patch to the local region
                inline bool isRegionPatch(const label patchI) const;

                //- Return the list of patch IDs on the primary region coupled
                //  to this region
                inline const labelList& primaryPatchIDs() const;

                //- Return the list of patch IDs internally coupled with the
                //  primary region
                inline const labelList& intCoupledPatchIDs() const;

                //- Return region ID corresponding to primaryPatchID
                inline label regionPatchID(const label primaryPatchID) const;


        // Helper

            //- Convert a local region field to the primary region
            template<class Type>
            void toPrimary
            (
                const label regionPatchI,
                List<Type>& regionField
            ) const;

            //- Convert a primary region field to the local region
            template<class Type>
            void toRegion
            (
                const label regionPatchI,
                List<Type>& primaryFieldField
            ) const;

            //- Convert a local region field to the primary region with op
            template<class Type, class BinaryOp>
            void toPrimary
            (
                const label regionPatchI,
                List<Type>& regionField,
                const BinaryOp& bop
            ) const;

            //- Convert a primary region field to the local region with op
            template<class Type, class BinaryOp>
            void toRegion
            (
                const label regionPatchI,
                List<Type>& primaryFieldField,
                const BinaryOp& bop
            ) const;


        // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the region
            virtual void evolveRegion();

            //- Evolve the film
            virtual void evolve();


        // I-O

            //- Provide some feedback
            virtual void info() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regionModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
