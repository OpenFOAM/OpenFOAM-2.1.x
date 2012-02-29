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
    Foam::interpolatedSolidThermo

Description
    Table interpolated solid thermo

SourceFiles
    interpolatedSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedSolidThermo_H
#define interpolatedSolidThermo_H

#include "basicSolidThermo.H"
#include "interpolateSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class interpolatedSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class interpolatedSolidThermo
:
    public basicSolidThermo,
    public interpolateSolid
{
private:

    // Private data

        //- Dictionary
        const dictionary dict_;


        // Helper functions

            //- Density [kg/m3]
            tmp<scalarField> rho(const label patchI) const;

            //- Scatter coefficient [1/m]
            tmp<scalarField> sigmaS(const label) const;

            //- Absorption coefficient [1/m]
            tmp<scalarField> kappa(const label) const;

             //- Emissivity []
            tmp<scalarField> emissivity(const label) const;


public:

    // Constructors

        //- Construct from mesh
        interpolatedSolidThermo(const fvMesh& mesh, const word);

        //- Construct from mesh and dictionary
        interpolatedSolidThermo
        (
            const fvMesh& mesh,
            const word,
            const dictionary& dict
        );



    //- Destructor
    virtual ~interpolatedSolidThermo();


    // Member Functions

        //- Calculate properties
        void calculate();


        // Derived properties

            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<volScalarField> Cp() const;

            //- Heat of formation [J/kg]
            virtual tmp<volScalarField> Hf() const;


        // Per patch calculation

            //- Specific heat capacity [J/kg/K)]
            virtual tmp<scalarField> Cp(const label patchI) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const;


        // I-O

            //- Write the interpolatedSolidThermo properties
            virtual bool writeData(Ostream& os) const;

            //- Read the interpolatedSolidThermo properties
            virtual bool read();

            //- Read the interpolatedSolidThermo properties
            bool read(const dictionary& dict);

            //- Ostream Operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const interpolatedSolidThermo& s
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
