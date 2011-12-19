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
    Foam::porousBafflePressureFvPatchField

Description
    Foam::porousBafflePressureFvPatchField
    the porous baffle operates on a cyclic patch and introduce a jump on the p
    field as follow:

    deltaP = -(I*mu*U + 0.5*D*rho*magSqr(U)*L)

    where:

    I is the inertial coefficient
    D is the darcy coeafficient
    L is the porous media lenght in the flow direction

SourceFiles
    porousBafflePressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef porousBafflePressureFvPatchField_H
#define porousBafflePressureFvPatchField_H

#include "fixedJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class porousBafflePressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class porousBafflePressureFvPatchField
:
    public fixedJumpFvPatchField<Type>
{
    // Private data

        //- Darcy pressure lost coefficient
        scalar D_;

        //- Inertia pressure lost coefficient
        scalar I_;

        //- Porous media length
        scalar length_;


public:

    //- Runtime type information
    TypeName("porousBafflePressure");


    // Constructors

        //- Construct from patch and internal field
        porousBafflePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousBafflePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given porousBafflePressureFvPatchField
        // onto a new patch
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new porousBafflePressureFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new porousBafflePressureFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


//- Specialisation of the jump-condition for the pressure
template<>
void porousBafflePressureFvPatchField<scalar>::updateCoeffs();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porousBafflePressureFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
