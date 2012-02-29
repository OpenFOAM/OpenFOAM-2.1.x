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

Class
    Foam::constSolidThermo

Description
    The thermophysical properties of a constSolidThermo

SourceFiles
    constSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constSolidThermo_H
#define constSolidThermo_H

#include "basicSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class constSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class constSolidThermo
:
    public basicSolidThermo
{
private:

    // Private data

        //- Dictionary
        dictionary dict_;

        //- Constant thermal conductivity [W/(m.K)]
        dimensionedScalar constK_;

        //- Thermal conductivity field[W/(m.K)]
        volScalarField K_;

        //- Density [kg/m3]
        dimensionedScalar constRho_;

        //- Specific heat capacity [J/(kg.K)]
        dimensionedScalar constCp_;

        //- Heat of formation [J/kg]
        dimensionedScalar constHf_;

        //- Emissivity
        dimensionedScalar constEmissivity_;

        //- Absorptivity [1/m]
        dimensionedScalar constKappa_;

        //- Scatter [1/m]
        dimensionedScalar constSigmaS_;


public:

    //- Runtime type information
    TypeName("constSolidThermo");


    // Constructors

        //- Construct from mesh
        constSolidThermo(const fvMesh& mesh);

         //- Construct from mesh and dict
        constSolidThermo(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~constSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


        // Derived thermal properties

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> K() const;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volSymmTensorField> directionalK() const;

            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<volScalarField> Cp() const;

            //- Heat of formation [J/kg]
            virtual tmp<volScalarField> Hf() const;


        // Per patch calculation

            //- Thermal conductivity [W//m/K]
            virtual tmp<scalarField> K(const label patchI) const;

            //- Thermal conductivity [W//m/K]
            virtual tmp<symmTensorField>directionalK(const label) const;

            //- Specific heat capacity [J/kg/K)]
            virtual tmp<scalarField> Cp(const label patchI) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const;



        // I-O

            //- Write the constSolidThermo properties
            virtual bool writeData(Ostream& os) const;

            //- Read solidThermophysicalProperties dictionary
            virtual bool read();

            //- Read solidThermophysicalProperties dictionary
            bool read(const dictionary&);

            //- Ostream Operator
            friend Ostream& operator<<(Ostream& os, const constSolidThermo& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
