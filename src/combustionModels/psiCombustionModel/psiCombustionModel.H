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
    Foam::psiCombustionModel

Description
    Combustion models for compressibility-based thermodynamics

SourceFiles
    psiCombustionModelI.H
    psiCombustionModel.C
    psiCombustionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef psiCombustionModel_H
#define psiCombustionModel_H

#include "combustionModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "hsCombustionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                     class psiCombustionModel Declaration
\*---------------------------------------------------------------------------*/

class psiCombustionModel
:
    public combustionModel
{

    // Private Member Functions

        //- Construct as copy (not implemented)
        psiCombustionModel(const psiCombustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const psiCombustionModel&);


protected:

    // Protected data


        //- Thermo package
        autoPtr<hsCombustionThermo> thermo_;




public:

    //- Runtime type information
    TypeName("psiCombustionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiCombustionModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh
        ),
        (modelType, mesh)
    );


    // Constructors


        //- Construct from components and thermo
        psiCombustionModel
        (
            const word& modelType,
            const fvMesh& mesh
        );



    //- Selector
    static autoPtr<psiCombustionModel> New
    (
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~psiCombustionModel();


    // Member Functions

        //- Return access to the thermo package
        inline hsCombustionThermo& thermo();

        //- Return const access to the thermo package
        inline const hsCombustionThermo& thermo() const;

        //- Return tmp of rho
        virtual tmp<volScalarField> rho() const;

    // I-O

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace combustionModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "psiCombustionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
