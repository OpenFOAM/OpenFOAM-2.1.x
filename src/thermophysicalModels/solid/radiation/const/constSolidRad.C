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

#include "constSolidRad.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
constSolidRad<thermo>::constSolidRad(const dictionary& dict)
:
    thermo(dict),
    kappa_(readScalar(dict.subDict("radiation").lookup("kappa"))),
    sigmaS_(readScalar(dict.subDict("radiation").lookup("sigmaS"))),
    emissivity_(readScalar(dict.subDict("radiation").lookup("emissivity")))
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class thermo>
Ostream& operator<<(Ostream& os, const constSolidRad<thermo>& pg)
{
    os  << static_cast<const thermo&>(pg);
    os << tab << pg.kappa_ << tab << pg.sigmaS_ << tab << pg.emissivity_;

    os.check("Ostream& operator<<(Ostream& os, const constSolidRad& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
