/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::incompressible::LESModels::nuSgsURoughWallFunctionFvPatchScalarField

Description
    Wall function boundary condition for rough walls, based on velocity.

SourceFiles
    nuSgsURoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nuSgsURoughWallFunctionFvPatchScalarField_H
#define nuSgsURoughWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
          Class nuSgsURoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nuSgsURoughWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Von Karman constant
        scalar kappa_;

        //- E constant
        scalar E_;

        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;

        // Roughness model parameters

            //- Height
            scalar roughnessHeight_;

            //- Constant
            scalar roughnessConstant_;

            //- Scale factor
            scalar roughnessFactor_;


    // Protected Member Functions

        //- Calculate the Y+ at the edge of the laminar sublayer
        scalar calcYPlusLam(const scalar kappa, const scalar E) const;

        //- Calculate yPLus
        tmp<scalarField> calcYPlus(const scalarField& magUp) const;


public:

    //- Runtime type information
    TypeName("nuSgsURoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nuSgsURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nuSgsURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nuSgsURoughWallFunctionFvPatchScalarField
        //  onto a new patch
        nuSgsURoughWallFunctionFvPatchScalarField
        (
            const nuSgsURoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nuSgsURoughWallFunctionFvPatchScalarField
        (
            const nuSgsURoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nuSgsURoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nuSgsURoughWallFunctionFvPatchScalarField
        (
            const nuSgsURoughWallFunctionFvPatchScalarField&,
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
                new nuSgsURoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the roughness height
            scalar roughnessHeight() const
            {
                return roughnessHeight_;
            }

            //- Return reference to the roughness height to allow adjustment
            scalar& roughnessHeight()
            {
                return roughnessHeight_;
            }


            //- Return the roughness constant scale
            scalar roughnessConstant() const
            {
                return roughnessConstant_;
            }

            //- Return reference to the roughness constant to allow adjustment
            scalar& roughnessConstant()
            {
                return roughnessConstant_;
            }

            //- Return the roughness scale factor
            scalar roughnessFactor() const
            {
                return roughnessFactor_;
            }

            //- Return reference to the roughness scale factor to allow
            //  adjustment
            scalar& roughnessFactor()
            {
                return roughnessFactor_;
            }


    // I-O

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
