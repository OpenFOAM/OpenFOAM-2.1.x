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
    Foam::psiChemistryCombustionModel

Description
    Combustion models for compressibility-based thermodynamics

SourceFiles
    psiChemistryCombustionModelI.H
    psiChemistryCombustionModel.C
    psiChemistryCombustionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef psiChemistryCombustionModel_H
#define psiChemistryCombustionModel_H

#include "combustionModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "psiChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                     class psiChemistryCombustionModel Declaration
\*---------------------------------------------------------------------------*/

class psiChemistryCombustionModel
:
    public combustionModel
{

    // Private Member Functions

        //- Construct as copy (not implemented)
        psiChemistryCombustionModel(const psiChemistryCombustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const psiChemistryCombustionModel&);


protected:

    // Protected data

        //- Auto pointer to psiChemistry
        autoPtr<psiChemistryModel> pChemistry_;


public:

    //- Runtime type information
    TypeName("psiChemistryCombustionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiChemistryCombustionModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh
        ),
        (modelType, mesh)
    );


    // Constructors


        //- Construct from components and thermo
        psiChemistryCombustionModel
        (
            const word& modelType,
            const fvMesh& mesh
        );



    //- Selector
    static autoPtr<psiChemistryCombustionModel> New
    (
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~psiChemistryCombustionModel();


    // Member Functions


        //- Return access to the thermo package
        inline psiChemistryModel& pChemistry();

        //- Return const access to the thermo package
        inline const psiChemistryModel& pChemistry() const;

        //- Return const access to rho
        inline tmp<volScalarField> rho() const;

        //- Return const access to rho
        inline const hsCombustionThermo& thermo() const;

        //- Return non const access to rho
        inline hsCombustionThermo& thermo();


    // I-O

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace combustionModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "psiChemistryCombustionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
