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
    Foam::phaseHydrostaticPressureFvPatchScalarField

Description
    Phase hydrostatic pressure boundary condition calculated as

    pRefValue + rho*g.(x - pRefPoint)

    where rho is provided and assumed uniform

    applied according to the phase-fraction field provided:
        1 -> fix value to that provided
        0 -> zero-gradient

SourceFiles
    phaseHydrostaticPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef phaseHydrostaticPressureFvPatchScalarField_H
#define phaseHydrostaticPressureFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class phaseHydrostaticPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseHydrostaticPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

        //- Name of phase-fraction field
        word phaseName_;

        //- Constant density in the far-field
        scalar rho_;

        //- Reference pressure
        scalar pRefValue_;

        //- Reference pressure location
        vector pRefPoint_;


public:

    //- Runtime type information
    TypeName("phaseHydrostaticPressure");


    // Constructors

        //- Construct from patch and internal field
        phaseHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        phaseHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  phaseHydrostaticPressureFvPatchScalarField onto a new patch
        phaseHydrostaticPressureFvPatchScalarField
        (
            const phaseHydrostaticPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        phaseHydrostaticPressureFvPatchScalarField
        (
            const phaseHydrostaticPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new phaseHydrostaticPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        phaseHydrostaticPressureFvPatchScalarField
        (
            const phaseHydrostaticPressureFvPatchScalarField&,
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
                new phaseHydrostaticPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the phaseName
            const word& phaseName() const
            {
                return phaseName_;
            }

            //- Return reference to the phaseName to allow adjustment
            word& phaseName()
            {
                return phaseName_;
            }

            //- Return the constant density in the far-field
            scalar rho() const
            {
                return rho_;
            }

            //- Return reference to the constant density in the far-field
            //  to allow adjustment
            scalar& rho()
            {
                return rho_;
            }

            //- Return the reference pressure
            scalar pRefValue() const
            {
                return pRefValue_;
            }

            //- Return reference to the reference pressure to allow adjustment
            scalar& pRefValue()
            {
                return pRefValue_;
            }

            //- Return the pressure reference location
            const vector& pRefPoint() const
            {
                return pRefPoint_;
            }

            //- Return reference to the pressure reference location
            //  to allow adjustment
            vector& pRefPoint()
            {
                return pRefPoint_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchScalarField& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
