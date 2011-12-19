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
    Foam::compressible::LESModels::
        muSgsUSpaldingWallFunctionFvPatchScalarField

Description
    Wall function boundary condition for walls, based on velocity, using
    Spaldings law to give a continuous muSgs profile to the wall (y+ = 0).

SourceFiles
    muSgsUSpaldingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef muSgsUSpaldingWallFunctionFvPatchScalarField_H
#define muSgsUSpaldingWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
       Class muSgsUSpaldingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class muSgsUSpaldingWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of velocity field
        word UName_;

        //- Name of density field
        word rhoName_;

        //- Name of laminar viscosity field
        word muName_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;


public:

    //- Runtime type information
    TypeName("muSgsUSpaldingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        muSgsUSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        muSgsUSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  muSgsUSpaldingWallFunctionFvPatchScalarField
        //  onto a new patch
        muSgsUSpaldingWallFunctionFvPatchScalarField
        (
            const muSgsUSpaldingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        muSgsUSpaldingWallFunctionFvPatchScalarField
        (
            const muSgsUSpaldingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new muSgsUSpaldingWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        muSgsUSpaldingWallFunctionFvPatchScalarField
        (
            const muSgsUSpaldingWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new muSgsUSpaldingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
