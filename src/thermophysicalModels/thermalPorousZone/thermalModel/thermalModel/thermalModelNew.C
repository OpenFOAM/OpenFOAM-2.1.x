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

#include "thermalModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::porousMedia::thermalModel>
Foam::porousMedia::thermalModel::New
(
    const porousZone& pZone
)
{
    // a missing thermalModel is the same as type "none"
    word modelType("none");

    if (const dictionary* dictPtr = pZone.dict().subDictPtr("thermalModel"))
    {
        dictPtr->lookup("type") >> modelType;
    }

    Info<< "Selecting thermalModel " << modelType << endl;

    pZoneConstructorTable::iterator cstrIter =
        pZoneConstructorTablePtr_->find(modelType);

    if (cstrIter == pZoneConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "porousMedia::thermalModel::New(const porousZone&)"
        )   << "Unknown thermalModel type "
            << modelType << nl << nl
            << "Valid thermalModel types are :" << endl
            << pZoneConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<porousMedia::thermalModel>(cstrIter()(pZone));
}


// ************************************************************************* //
