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
    Foam::thermoBaffle2D

Description
    2D thermal baffle

SourceFiles
    thermoBaffle2D.C
    thermoBaffle2DI.H

\*---------------------------------------------------------------------------*/

#ifndef thermoBaffle2D_H
#define thermoBaffle2D_H

#include "thermoBaffleModel.H"
#include "volFieldsFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermoBaffleModels
{


/*---------------------------------------------------------------------------*\
                        Class thermoBaffle2D Declaration
\*---------------------------------------------------------------------------*/

class thermoBaffle2D
:
    public thermoBaffleModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoBaffle2D(const thermoBaffle2D&);

        //- Disallow default bitwise assignment
        void operator=(const thermoBaffle2D&);

        //- Initialize thermoBaffle2D
        void init();


protected:

    // Protected data

        // Solution parameters

            //- Number of non orthogonal correctors
            label nNonOrthCorr_;


        // Thermo properties

            //- Solid thermo
            autoPtr<basicSolidThermo> thermo_;

            //- Temperature / [K]
            volScalarField& T_;


        // Source term fields

            //- Surface energy source  / [J/m2/s]
            volScalarField Qs_;

            //- Volumetric energy source  / [J/m3/s]
            volScalarField Q_;


    // Protected member functions

        //- Read control parameters from IO dictionary
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);


        // Equations

            //- Solve energy equation
            void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermoBaffle2D");


    // Constructors

        //- Construct from components
        thermoBaffle2D(const word& modelType, const fvMesh& mesh);

        //- Construct from components and dict
        thermoBaffle2D
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict

        );



    //- Destructor
    virtual ~thermoBaffle2D();


    // Member Functions

        // Thermo properties

            //- Return const reference to the basicSolidThermo
            virtual const basicSolidThermo& thermo() const;


            // Fields

                //- Return the film specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const;

                //- Return solid absortivity [1/m]
                virtual const volScalarField& kappa() const;

                //- Return the film mean temperature [K]
                virtual const volScalarField& T() const;

                //- Return density [Kg/m3]
                virtual const volScalarField& rho() const;

                //- Return thermal conductivity [W/m/K]
                virtual const volScalarField& K() const;


            // Helper functions

                //- Return sensible enthalpy as a function of temperature
                //  for a patch
                inline tmp<scalarField> hs
                (
                    const scalarField& T,
                    const label patchI
                ) const;

                //- Return sensible enthalpy as a function of temperature
                inline tmp<volScalarField> hs
                (
                    const volScalarField& T
                ) const;

                //- Return temperature as a function of sensible enthalpy
                inline tmp<volScalarField> T
                (
                    const volScalarField& hs
                ) const;


        // Evolution

            //- Pre-evolve  thermal baffle
            virtual void preEvolveRegion();

            //- Evolve the thermal baffle
            virtual void evolveRegion();



       // I-O

            //- Provide some feedback
            virtual void info() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoBaffleModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoBaffle2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
