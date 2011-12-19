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
    Foam::compressible::RASModels::mutkRoughWallFunctionFvPatchScalarField

Description
    Boundary condition for turbulent (kinematic) viscosity when using wall
    functions for rough walls, based on turbulence kinetic energy.

    Manipulates the E parameter to account for roughness effects, based on
    KsPlus.

    - roughness height = sand-grain roughness (0 for smooth walls)
    - roughness constant = 0.5-1.0 (0.5 default)

SourceFiles
    mutkRoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutkRoughWallFunctionFvPatchScalarField_H
#define compressibleMutkRoughWallFunctionFvPatchScalarField_H

#include "mutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
          Class mutkRoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutkRoughWallFunctionFvPatchScalarField
:
    public mutkWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- Roughness height
        scalarField Ks_;

        //- Roughness constant
        scalarField Cs_;


    // Private Member Functions

        //- Compute the roughness function
        virtual scalar fnRough(const scalar KsPlus, const scalar Cs) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcMut() const;


public:

    //- Runtime type information
    TypeName("mutkRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutkRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutkRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutkRoughWallFunctionFvPatchScalarField
        //  onto a new patch
        mutkRoughWallFunctionFvPatchScalarField
        (
            const mutkRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutkRoughWallFunctionFvPatchScalarField
        (
            const mutkRoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mutkRoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mutkRoughWallFunctionFvPatchScalarField
        (
            const mutkRoughWallFunctionFvPatchScalarField&,
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
                new mutkRoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
