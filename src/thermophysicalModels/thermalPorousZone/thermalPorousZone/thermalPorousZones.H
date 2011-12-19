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

Class
    Foam::thermalPorousZones

Description
    A centralized thermalPorousZone collection.

    Container class for a set of thermalPorousZones with the thermalPorousZone
    member functions implemented to loop over the functions for each
    thermalPorousZone.

    The input file \c constant/thermalPorousZone is implemented as
    IOPtrList\<thermalPorousZone\> (but written as a dictionary)
    and contains the following type of data:

    \verbatim
    1
    (
    cat1
    {
        coordinateSystem    system_10;
        porosity    0.781;
        Darcy
        {
            d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
            f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
        }
        thermalModel
        {
            type    fixedTemperature;
            T       T [0 0 1 0 0] 600;
        }
    }
    )
    \endverbatim

SourceFiles
    thermalPorousZones.C

\*---------------------------------------------------------------------------*/

#ifndef thermalPorousZones_H
#define thermalPorousZones_H

#include "PorousZones.H"
#include "thermalPorousZone.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class thermalPorousZones Declaration
\*---------------------------------------------------------------------------*/

class thermalPorousZones
:
    public PorousZones<thermalPorousZone>
{

public:

    // Constructors

        //- Construct from fvMesh
        thermalPorousZones(const fvMesh&);


    // Member Functions

        //- Add the thermal source to the enthalpy equation
        void addEnthalpySource
        (
            const basicThermo&,
            const volScalarField& rho,
            fvScalarMatrix& hEqn
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
