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

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DevolatilisationModel<CloudType>::DevolatilisationModel
(
    CloudType& owner
)
:
    SubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::DevolatilisationModel<CloudType>::DevolatilisationModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, type)
{}


template<class CloudType>
Foam::DevolatilisationModel<CloudType>::DevolatilisationModel
(
    const DevolatilisationModel<CloudType>& dm
)
:
    SubModelBase<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DevolatilisationModel<CloudType>::~DevolatilisationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DevolatilisationModel<CloudType>::calculate
(
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalarField&,
    bool&,
    scalarField&
) const
{
    notImplemented
    (
        "void Foam::DevolatilisationModel<CloudType>::calculate"
        "("
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalarField&, "
            "bool&, "
            "scalarField&"
        ") const"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DevolatilisationModelNew.C"

// ************************************************************************* //
