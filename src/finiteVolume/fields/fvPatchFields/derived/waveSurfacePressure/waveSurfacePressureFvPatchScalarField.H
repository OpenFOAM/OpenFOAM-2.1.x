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
    Foam::waveSurfacePressureFvPatchScalarField

Description
    Applies the surface wave pressure, based on the wave height

        p = -rho.g.zeta

    Where

        p = pressure - kinematic of dynamic depending on the flux units
        zeta = wave height vector [m]
        g = acceleration due to gravity [m/s2]


    Note:

        This boundary also updates the wave height boundary field, which
        is accessed via lookup from the database


SourceFiles
    waveSurfacePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef waveSurfacePressureFvPatchScalarField_H
#define waveSurfacePressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class waveSurfacePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waveSurfacePressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Flux field name
        word phiName_;

        //- Density field name (if compressible)
        word rhoName_;

        //- Wave height field name
        word zetaName_;

        //- Old-time zeta patch value
        vectorField zeta0_;

        //- Current time index used to store ms0_
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("waveSurfacePressure");


    // Constructors

        //- Construct from patch and internal field
        waveSurfacePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveSurfacePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given waveSurfacePressureFvPatchScalarField
        //  onto a new patch
        waveSurfacePressureFvPatchScalarField
        (
            const waveSurfacePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waveSurfacePressureFvPatchScalarField
        (
            const waveSurfacePressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new waveSurfacePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveSurfacePressureFvPatchScalarField
        (
            const waveSurfacePressureFvPatchScalarField&,
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
                new waveSurfacePressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
