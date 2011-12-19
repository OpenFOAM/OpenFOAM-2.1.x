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
    Foam::solidChemistryModel

Description
    Chemistry model for solid thermodynamics

SourceFiles
    solidChemistryModelI.H
    solidChemistryModel.C
    newChemistrySolidModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidChemistryModel_H
#define solidChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "basicSolidThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class solidChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class solidChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        solidChemistryModel(const solidChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidChemistryModel&);


protected:

    // Protected data

        //- Solid thermo package
        autoPtr<basicSolidThermo> solidThermo_;


public:

    //- Runtime type information
    TypeName("solidChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& solidThermoTypeName
        ),
        (mesh, compTypeName, solidThermoTypeName)
    );


    // Constructors

        //- Construct from mesh and thermo type name
        solidChemistryModel
        (
            const fvMesh& mesh,
            const word& solidThermoTypeName
        );


    //- Selector
    static autoPtr<solidChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~solidChemistryModel();


    // Member Functions

        //- Return access to the solid thermo package
        inline basicSolidThermo& solidThermo();

        //- Return const access to the solid thermo package
        inline const basicSolidThermo& solidThermo() const;

        //- Return total gases mass source term [kg/m3/s]
        virtual tmp<volScalarField> RRg() const = 0;

        //- Return total solids mass source term [kg/m3/s]
        virtual tmp<volScalarField> RRs() const = 0;

        //- Return chemical source terms for solids [kg/m3/s]
        virtual tmp<volScalarField> RRs(const label i) const = 0;

        //- Return chemical source terms for gases [kg/m3/s]
        virtual tmp<volScalarField> RRg(const label i) const = 0;

        //- Return sensible enthalpy for gas i [J/Kg]
        virtual tmp<volScalarField> gasHs
        (
            const volScalarField& T,
            const label i
        ) const = 0;

        //- Return specie Table for gases
        virtual const speciesTable& gasTable() const = 0;

        //- Set reacting status of cell, cellI
        virtual void setCellReacting(const label cellI, const bool active) = 0;

        //- Calculates the reaction rates
        virtual void calculate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
