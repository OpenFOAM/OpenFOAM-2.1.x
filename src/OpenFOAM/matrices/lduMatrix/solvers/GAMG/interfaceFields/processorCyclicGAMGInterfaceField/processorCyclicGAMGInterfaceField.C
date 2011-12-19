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

#include "processorCyclicGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorCyclicGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        processorCyclicGAMGInterfaceField,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorCyclicGAMGInterfaceField::processorCyclicGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    processorGAMGInterfaceField(GAMGCp, fineInterface)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorCyclicGAMGInterfaceField::~processorCyclicGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//void Foam::processorCyclicGAMGInterfaceField::initInterfaceMatrixUpdate
//(
//    const scalarField& psiInternal,
//    scalarField&,
//    const lduMatrix&,
//    const scalarField&,
//    const direction,
//    const Pstream::commsTypes commsType
//) const
//{
//    procInterface_.compressedSend
//    (
//        commsType,
//        procInterface_.interfaceInternalField(psiInternal)()
//    );
//}
//
//
//void Foam::processorCyclicGAMGInterfaceField::updateInterfaceMatrix
//(
//    const scalarField&,
//    scalarField& result,
//    const lduMatrix&,
//    const scalarField& coeffs,
//    const direction cmpt,
//    const Pstream::commsTypes commsType
//) const
//{
//    scalarField pnf
//    (
//        procInterface_.compressedReceive<scalar>(commsType, coeffs.size())
//    );
//    transformCoupleField(pnf, cmpt);
//
//    const labelUList& faceCells = procInterface_.faceCells();
//
//    forAll(faceCells, elemI)
//    {
//        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
//    }
//}


// ************************************************************************* //
