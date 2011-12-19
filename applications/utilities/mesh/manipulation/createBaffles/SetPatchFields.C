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

#include "SetPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class GeoField>
void Foam::SetPatchFields
(
    PtrList<GeoField>& fields,
    const label patchI,
    const typename GeoField::value_type& initVal
)
{
    forAll(fields, i)
    {
        typename GeoField::PatchFieldType& pfld =
            fields[i].boundaryField()[patchI];
        pfld == initVal;
    }
}


// ************************************************************************* //
