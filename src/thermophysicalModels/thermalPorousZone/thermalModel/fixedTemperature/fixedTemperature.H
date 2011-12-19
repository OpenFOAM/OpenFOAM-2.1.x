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
    Foam::porousMedia::fixedTemperature

Description
    Fixed temperature model for porous media

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperature_H
#define fixedTemperature_H

#include "thermalModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousMedia
{

/*---------------------------------------------------------------------------*\
                      Class fixedTemperature Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperature
:
    public thermalModel
{

protected:

    // Protected data

        //- Fixed temperature
        const scalar T_;


public:

    //- Runtime type information
    TypeName("fixedTemperature");


    // Constructors

        //- Construct from porous zone
        fixedTemperature(const porousZone& pZone);


    //- Destructor
    virtual ~fixedTemperature();


    // Member Functions

        //- Add the thermal source to the enthalpy equation
        virtual void addEnthalpySource
        (
            const basicThermo&,
            const volScalarField& rho,
            fvScalarMatrix& hEqn
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousMedia

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

