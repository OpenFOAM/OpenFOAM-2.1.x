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
    Foam::thermoBaffleModel

Description

SourceFiles
    thermoBaffleModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermoBaffleModel_H
#define thermoBaffleModel_H

#include "runTimeSelectionTables.H"
#include "scalarIOField.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "basicSolidThermo.H"
#include "regionModel1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermoBaffleModels
{

/*---------------------------------------------------------------------------*\
                      Class thermoBaffleModel Declaration
\*---------------------------------------------------------------------------*/

class thermoBaffleModel
:
    public regionModel1D
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermoBaffleModel(const thermoBaffleModel&);

        //- Disallow default bitwise assignment
        void operator=(const thermoBaffleModel&);

        //- Initialize thermal Baffle
        void init();


protected:

    // Protected Data

        //- Baffle physical thickness
        scalarField thickness_;

        //- Baffle mesh thickness
        dimensionedScalar delta_;

        //- Is it one dimension
        bool oneD_;

        //- Is thickness constant
        bool constantThickness_;


    // Protected Member Functions

        //- Read control parameters from IOdictionary
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary&);


public:

    //- Runtime type information
    TypeName("thermoBaffleModel");


    // Declare runtime constructor selection tables

         declareRunTimeSelectionTable
         (
             autoPtr,
             thermoBaffleModel,
             mesh,
             (
                const word& modelType,
                const fvMesh& mesh
             ),
             (modelType, mesh)
         );

         declareRunTimeSelectionTable
         (
             autoPtr,
             thermoBaffleModel,
             dictionary,
             (
                const word& modelType,
                const fvMesh& mesh,
                const dictionary& dict
             ),
             (modelType, mesh, dict)
         );


    // Constructors

        //- Construct null from mesh
        thermoBaffleModel(const fvMesh& mesh);

        //- Construct from type name and mesh
        thermoBaffleModel(const word& modelType, const fvMesh& mesh);

        //- Construct from type name and mesh and dict
        thermoBaffleModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected model
         static autoPtr<thermoBaffleModel> New(const fvMesh& mesh);

        //- Return a reference to the selected model using dictionary
        static autoPtr<thermoBaffleModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermoBaffleModel();


    // Member Functions


        // Access

            //- Return solid thermo
            virtual const basicSolidThermo& thermo() const = 0;

            //- Return thickness
            const scalarField& thickness() const
            {
                return thickness_;
            }

            //- Return geometrical thickness
            const dimensionedScalar& delta() const
            {
                return delta_;
            }

             //- Return if region is one dimensional
            bool oneD() const
            {
                return oneD_;
            }

             //- Return if region has constant thickness
            bool constantThickness() const
            {
                return constantThickness_;
            }


            // Fields

                //- Return density [kg/m3]
                virtual const volScalarField& rho() const = 0;

                //- Return const temperature [K]
                virtual const volScalarField& T() const = 0;

                //- Return specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const = 0;

                //- Return the region absorptivity [1/m]
                virtual const volScalarField& kappa() const = 0;

                //- Return the region thermal conductivity [W/m/k]
                virtual const volScalarField& K() const = 0;


        // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoBaffleModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
