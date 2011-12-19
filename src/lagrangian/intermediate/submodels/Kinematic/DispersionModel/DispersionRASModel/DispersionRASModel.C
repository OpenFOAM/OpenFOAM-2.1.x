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

#include "DispersionRASModel.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DispersionRASModel<CloudType>::DispersionRASModel
(
    const dictionary&,
    CloudType& owner
)
:
    DispersionModel<CloudType>(owner),
    turbulence_
    (
        owner.mesh().objectRegistry::template lookupObject
        <
            compressible::RASModel
        >
        (
            "RASProperties"
        )
    ),
    kPtr_(NULL),
    ownK_(false),
    epsilonPtr_(NULL),
    ownEpsilon_(false)
{}


template<class CloudType>
Foam::DispersionRASModel<CloudType>::DispersionRASModel
(
    DispersionRASModel<CloudType>& dm
)
:
    DispersionModel<CloudType>(dm),
    turbulence_(dm.turbulence_),
    kPtr_(dm.kPtr_),
    ownK_(dm.ownK_),
    epsilonPtr_(dm.epsilonPtr_),
    ownEpsilon_(dm.ownEpsilon_)
{
    dm.ownK_ = false;
    dm.ownEpsilon_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DispersionRASModel<CloudType>::~DispersionRASModel()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DispersionRASModel<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        tmp<volScalarField> tk = this->turbulence().k();
        if (tk.isTmp())
        {
            kPtr_ = tk.ptr();
            ownK_ = true;
        }
        else
        {
            kPtr_ = tk.operator->();
            ownK_ = false;
        }

        tmp<volScalarField> tepsilon = this->turbulence().epsilon();
        if (tepsilon.isTmp())
        {
            epsilonPtr_ = tepsilon.ptr();
            ownEpsilon_ = true;
        }
        else
        {
            epsilonPtr_ = tepsilon.operator->();
            ownEpsilon_ = false;
        }
    }
    else
    {
        if (ownK_ && kPtr_)
        {
            deleteDemandDrivenData(kPtr_);
            ownK_ = false;
        }
        if (ownEpsilon_ && epsilonPtr_)
        {
            deleteDemandDrivenData(epsilonPtr_);
            ownEpsilon_ = false;
        }
    }
}


template<class CloudType>
void Foam::DispersionRASModel<CloudType>::write(Ostream& os) const
{
    DispersionModel<CloudType>::write(os);

    os.writeKeyword("ownK") << ownK_ << token::END_STATEMENT << endl;
    os.writeKeyword("ownEpsilon") << ownEpsilon_ << token::END_STATEMENT
        << endl;
}


// ************************************************************************* //
