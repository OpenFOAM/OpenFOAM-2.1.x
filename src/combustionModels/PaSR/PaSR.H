/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    Foam::combustionModels::PaSR

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    PaSR.C

\*---------------------------------------------------------------------------*/

#ifndef PaSR_H
#define PaSR_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class PaSR Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType>
class PaSR
:
    public CombThermoType
{
    // Private data

        //- Mixing constant
        dimensionedScalar Cmix_;

        //- Turbulent reaction switch
        Switch turbulentReaction_;

        //- Mixing parameter
        volScalarField kappa_;

        //- Use reaction rate
        bool useReactionRate_;


    // Private Member Functions

        //- Return the chemical time scale
        tmp<volScalarField> tc() const;

        //- Disallow copy construct
        PaSR(const PaSR&);

        //- Disallow default bitwise assignment
        void operator=(const PaSR&);


public:

    //- Runtime type information
    TypeName("PaSR");


    // Constructors

        //- Construct from components
        PaSR
        (
            const word& modelType,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~PaSR();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix.
            virtual tmp<fvScalarMatrix> R(const volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //-  Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

    // I-O

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PaSR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
