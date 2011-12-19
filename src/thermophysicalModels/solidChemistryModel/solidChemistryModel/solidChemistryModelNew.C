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

#include "solidChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidChemistryModel> Foam::solidChemistryModel::New
(
    const fvMesh& mesh
)
{
    IOdictionary chemistryPropertiesDict
    (
        IOobject
        (
            "chemistryProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word userModel(chemistryPropertiesDict.lookup("solidChemistryModel"));

    const word ODEModelName(chemistryPropertiesDict.lookup("chemistrySolver"));
    const word gasThermoName(chemistryPropertiesDict.lookup("gasThermoModel"));

    // construct chemistry model type name by inserting first template argument
    const label tempOpen = userModel.find('<');
    const label tempClose = userModel.find('>');

    const word className = userModel(0, tempOpen);
    const word thermoTypeName =
        userModel(tempOpen + 1, tempClose - tempOpen - 1);

    const word modelType =
        ODEModelName + '<' + className
      + '<' + typeName + ',' + thermoTypeName + ',' + gasThermoName + ">>";


    if (debug)
    {
        Info<< "Selecting solidChemistryModel " << modelType << endl;
    }
    else
    {
        Info<< "Selecting solidChemistryModel " << userModel + gasThermoName
            << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(modelType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        if (debug)
        {
            FatalErrorIn("solidChemistryModel::New(const mesh&)")
                << "Unknown solidChemistryModel type "
                << modelType << nl << nl
                << "Valid solidChemistryModel types are:" << nl
                << fvMeshConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }
        else
        {
            wordList models = fvMeshConstructorTablePtr_->sortedToc();
            forAll(models, i)
            {
                models[i] = models[i].replace(typeName + ',', "");
            }

            FatalErrorIn("solidChemistryModel::New(const mesh&)")
                << "Unknown solidChemistryModel type "
                << userModel << nl << nl
                << "Valid solidChemistryModel types are:" << nl
                << models << nl
                << exit(FatalError);
        }
    }

    return autoPtr<solidChemistryModel>
        (cstrIter()(mesh, ODEModelName, thermoTypeName));
}


// ************************************************************************* //
