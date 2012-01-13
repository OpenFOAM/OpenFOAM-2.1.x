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
    baseName_("none"),
    name_("none"),
    coeffDict_(dictionary::null)
{}


template<class CloudType>
Foam::SubModelBase<CloudType>::SubModelBase
(
    CloudType& owner,
    const dictionary& dict,
    const word& baseName,
    const word& name,
    const word& dictExt
)
:
    owner_(owner),
    dict_(dict),
    baseName_(baseName),
    name_(name),
    coeffDict_(dict.subDict(name + dictExt))
{}


template<class CloudType>
Foam::SubModelBase<CloudType>::SubModelBase(const SubModelBase<CloudType>& smb)
:
    owner_(smb.owner_),
    dict_(smb.dict_),
    baseName_(smb.baseName_),
    name_(smb.name_),
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
const Foam::word& Foam::SubModelBase<CloudType>::baseName() const
{
    return baseName_;
}


template<class CloudType>
const Foam::word& Foam::SubModelBase<CloudType>::name() const
{
    return name_;
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
template<class Type>
Type Foam::SubModelBase<CloudType>::getBaseProperty
(
    const word& entryName
) const
{
    Type result = pTraits<Type>::zero;

    const dictionary& properties = this->owner().outputProperties();

    if (properties.found(baseName_))
    {
        const dictionary& baseDict = properties.subDict(baseName_);
        baseDict.readIfPresent(entryName, result);
    }

    return result;
}


template<class CloudType>
template<class Type>
void Foam::SubModelBase<CloudType>::setBaseProperty
(
    const word& entryName,
    const Type& value
)
{
    dictionary& properties = this->owner().outputProperties();

    if (properties.found(baseName_))
    {
        dictionary& baseDict = properties.subDict(baseName_);
        baseDict.add(entryName, value, true);
    }
    else
    {
        properties.add(baseName_, dictionary());
        properties.subDict(baseName_).add(entryName, value);
    }
}


template<class CloudType>
template<class Type>
Type Foam::SubModelBase<CloudType>::getModelProperty
(
    const word& entryName
) const
{
    Type result = pTraits<Type>::zero;

    const dictionary& properties = this->owner().outputProperties();

    if (properties.found(baseName_))
    {
        const dictionary& baseDict = properties.subDict(baseName_);

        if (baseDict.found(name_))
        {
            baseDict.subDict(name_).readIfPresent(entryName, result);
        }
    }

    return result;
}


template<class CloudType>
template<class Type>
void Foam::SubModelBase<CloudType>::getModelProperty
(
    const word& entryName,
    Type& value
) const
{
    const dictionary& properties = this->owner().outputProperties();

    if (properties.found(baseName_))
    {
        const dictionary& baseDict = properties.subDict(baseName_);

        if (baseDict.found(name_))
        {
            baseDict.subDict(name_).readIfPresent(entryName, value);
        }
    }
}


template<class CloudType>
template<class Type>
void Foam::SubModelBase<CloudType>::setModelProperty
(
    const word& entryName,
    const Type& value
)
{
    dictionary& properties = this->owner().outputProperties();

    if (properties.found(baseName_))
    {
        dictionary& baseDict = properties.subDict(baseName_);
        if (baseDict.found(name_))
        {
            baseDict.subDict(name_).add(entryName, value, true);
        }
        else
        {
            baseDict.add(name_, dictionary());
            baseDict.subDict(name_).add(entryName, value, true);
        }
    }
    else
    {
        properties.add(baseName_, dictionary());
        properties.subDict(baseName_).add(name_, dictionary());
        properties.subDict(baseName_).subDict(name_).add(entryName, value);
    }
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
