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

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicFvPatchScalarField_H
#define processorCyclicFvPatchScalarField_H

#include "processorCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<>
//void processorCyclicFvPatchField<scalar>::initInterfaceMatrixUpdate
//(
//    const scalarField&,
//    scalarField&,
//    const lduMatrix&,
//    const scalarField&,
//    const direction,
//    const Pstream::commsTypes commsType
//) const;
//
//
//template<>
//void processorCyclicFvPatchField<scalar>::updateInterfaceMatrix
//(
//    const scalarField&,
//    scalarField& result,
//    const lduMatrix&,
//    const scalarField& coeffs,
//    const direction,
//    const Pstream::commsTypes commsType
//) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
