/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "pressurePerfectGas.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressurePerfectGas::pressurePerfectGas(Istream& is)
:
    specie(is),
    pRef_(readScalar(is))
{
    is.check("pressurePerfectGas::pressurePerfectGas(Istream& is)");
}


Foam::pressurePerfectGas::pressurePerfectGas(const dictionary& dict)
:
    specie(dict),
    pRef_(readScalar(dict.subDict("equationOfState").lookup("pRef")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressurePerfectGas::write(Ostream& os) const
{
    specie::write(os);
    dictionary dict("equationOfState");
    dict.add("pRef", pRef_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pressurePerfectGas& pg)
{
    os  << static_cast<const specie&>(pg)
        << token::SPACE << pg.pRef_;

    os.check("Ostream& operator<<(Ostream& os, const pressurePerfectGas& st)");
    return os;
}


// ************************************************************************* //
