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

#include "specieThermo.H"
#include "IOstreams.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

template<class Thermo>
const Foam::scalar Foam::specieThermo<Thermo>::tol_ = 1.0e-4;

template<class Thermo>
const int Foam::specieThermo<Thermo>::maxIter_ = 100;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::specieThermo<Thermo>::specieThermo(Istream& is)
:
    Thermo(is)
{
    is.check("specieThermo<Thermo>::specieThermo(Istream&)");
}


template<class Thermo>
Foam::specieThermo<Thermo>::specieThermo(const dictionary& dict)
:
    Thermo(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::specieThermo<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const specieThermo<Thermo>& st)
{
    os  << static_cast<const Thermo&>(st);

    os.check("Ostream& operator<<(Ostream&, const specieThermo&)");
    return os;
}


// ************************************************************************* //
