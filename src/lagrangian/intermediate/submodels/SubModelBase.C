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

#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SubModelBase<CloudType>::SubModelBase(CloudType& owner)
:
    owner_(owner),
    dict_(dictionary::null),
    coeffDict_(dictionary::null)
{}


template<class CloudType>
Foam::SubModelBase<CloudType>::SubModelBase
(
    CloudType& owner,
    const dictionary& dict,
    const word& name,
    const word& dictExt
)
:
    owner_(owner),
    dict_(dict),
    coeffDict_(dict.subDict(name + dictExt))
{}


template<class CloudType>
Foam::SubModelBase<CloudType>::SubModelBase(const SubModelBase<CloudType>& smb)
:
    owner_(smb.owner_),
    dict_(smb.dict_),
    coeffDict_(smb.coeffDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SubModelBase<CloudType>::~SubModelBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::SubModelBase<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::SubModelBase<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::SubModelBase<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
bool Foam::SubModelBase<CloudType>::defaultCoeffs(const bool printMsg) const
{
    bool def = coeffDict_.lookupOrDefault<bool>("defaultCoeffs", false);
    if (printMsg && def)
    {
        Info<< incrIndent;
        Info<< indent << "Employing default coefficients" << endl;
        Info<< decrIndent;
    }

    return def;
}


template<class CloudType>
CloudType& Foam::SubModelBase<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
bool Foam::SubModelBase<CloudType>::active() const
{
    return true;
}


template<class CloudType>
void Foam::SubModelBase<CloudType>::cacheFields(const bool)
{
    // do nothing
}


template<class CloudType>
void Foam::SubModelBase<CloudType>::write(Ostream& os) const
{
    os.writeKeyword("owner") << owner_.name() << token::END_STATEMENT << nl;

    // not writing complete cloud dictionary, only coeffs
//    os  << dict_;
    os  << coeffDict_;
}


// ************************************************************************* //
