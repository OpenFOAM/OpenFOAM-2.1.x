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

#include "basicSourceList.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSourceList, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::basicSourceList::checkApplied() const
{
    if (mesh_.time().timeIndex() == checkTimeIndex_)
    {
        forAll(*this, i)
        {
            const basicSource& bs = this->operator[](i);
            bs.checkApplied();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSourceList::basicSourceList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<basicSource>(),
    mesh_(mesh),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{
    reset(dict);
}


Foam::basicSourceList::basicSourceList(const fvMesh& mesh)
:
    PtrList<basicSource>(),
    mesh_(mesh),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicSourceList::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        // safety:
        if (iter().isDict())
        {
            count ++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& sourceDict = iter().dict();

            this->set
            (
                i++,
                basicSource::New(name, sourceDict, mesh_)
            );
        }
    }
}


bool Foam::basicSourceList::read(const dictionary& dict)
{
    checkTimeIndex_ = mesh_.time().timeIndex() + 2;

    bool allOk = true;
    forAll(*this, i)
    {
        basicSource& bs = this->operator[](i);
        bool ok = bs.read(dict.subDict(bs.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::basicSourceList::writeData(Ostream& os) const
{
    // Write list contents
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    // Check state of IOstream
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const basicSourceList& sources
)
{
    sources.writeData(os);
    return os;
}


// ************************************************************************* //
