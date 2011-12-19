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
    Foam::basicSolidMixture

Description
    Foam::basicSolidMixture

\*---------------------------------------------------------------------------*/

#ifndef basicSolidMixture_H
#define basicSolidMixture_H

#include "volFields.H"
#include "speciesTable.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicSolidMixture Declaration
\*---------------------------------------------------------------------------*/

class basicSolidMixture
{

protected:

    // Protected data

        typedef speciesTable solidsTable;

        //- The names of the solids
        solidsTable components_;

        //- Solid mass fractions
        PtrList<volScalarField> Y_;


public:

    // Constructors

        //- Construct from word list and mesh
        basicSolidMixture
        (
            const wordList& solidNames,
            const fvMesh&
        );


    //- Destructor
    virtual ~basicSolidMixture()
    {}


    // Member Functions

        //- Return the solid table
        const solidsTable& components() const
        {
            return components_;
        }

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;


    // Derived cell based properties.

        //- Density
        virtual scalar rho(scalar T, label celli) const = 0;

        //- Absorption coefficient
        virtual scalar kappa(scalar T, label celli) const = 0;

        //- Scatter coefficient
        virtual scalar sigmaS(scalar T, label celli) const = 0;

        //- Thermal conductivity
        virtual scalar K(scalar T, label celli) const = 0;

        //- Emissivity coefficient
        virtual scalar emissivity(scalar T, label celli) const = 0;

        //- Formation enthalpy
        virtual scalar hf(scalar T, label celli) const = 0;

        //- Sensible enthalpy
        virtual scalar hs(scalar T, label celli) const = 0;

        //- Total enthalpy
        virtual scalar h(scalar T, label celli) const = 0;

        //- Specific heat capacity
        virtual scalar Cp(scalar T, label celli) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "basicSolidMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
