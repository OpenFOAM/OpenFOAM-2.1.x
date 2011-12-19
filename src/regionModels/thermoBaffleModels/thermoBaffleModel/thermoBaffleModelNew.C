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

#include "thermoBaffleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermoBaffleModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<thermoBaffleModel> thermoBaffleModel::New(const fvMesh& mesh)
{

    word modelType;
    {
        IOdictionary thermoBafflePropertiesDict
        (
            IOobject
            (
                "thermoBaffleProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        thermoBafflePropertiesDict.lookup("thermoBaffleModel") >> modelType;
    }
    Info<< "Selecting baffle model " << modelType << endl;

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(modelType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {

        FatalErrorIn("thermoBaffleModel::New(const fvMesh&)")
            << "Unknown thermoBaffleModel type " << modelType
            << nl << nl
            <<  "Valid thermoBaffleModel types are:" << nl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermoBaffleModel>(cstrIter()(modelType, mesh));
}


autoPtr<thermoBaffleModel> thermoBaffleModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{

    word modelType = dict.lookup("thermoBaffleModel");

    Info<< "Selecting baffle model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {

        FatalErrorIn("thermoBaffleModel::New(const fvMesh&,const dictionary&)")
            << "Unknown thermoBaffleModel type " << modelType
            << nl << nl
            <<  "Valid thermoBaffleModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermoBaffleModel>(cstrIter()(modelType, mesh, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoBaffleModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
