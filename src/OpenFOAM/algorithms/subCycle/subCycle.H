/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::subCycle

Description
    Perform a subCycleTime on a field

\*---------------------------------------------------------------------------*/

#ifndef subCycle_H
#define subCycle_H

#include "subCycleTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class subCycleField Declaration
\*---------------------------------------------------------------------------*/

template<class GeometricField>
class subCycleField
{
    // Private data

        //- Reference to the field being sub-cycled
        GeometricField& gf_;

        //- Copy of the "real" old-time value of the field
        GeometricField  gf0_;


public:

    // Constructors

        //- Construct field and number of sub-cycles
        subCycleField(GeometricField& gf)
        :
            gf_(gf),
            gf0_(gf.oldTime())
        {}


    //- Destructor
    ~subCycleField()
    {
        // Correct the time index of the field to correspond to the global time
        gf_.timeIndex() = gf_.time().timeIndex();

        // Reset the old-time field value
        gf_.oldTime() = gf0_;
        gf_.oldTime().timeIndex() = gf0_.timeIndex();
    }
};


/*---------------------------------------------------------------------------*\
                          Class subCycle Declaration
\*---------------------------------------------------------------------------*/

template<class GeometricField>
class subCycle
:
    public subCycleField<GeometricField>,
    public subCycleTime
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        subCycle(const subCycle<GeometricField>&);

        //- Disallow default bitwise assignment
        void operator=(const subCycle<GeometricField>&);


public:

    // Constructors

        //- Construct field and number of sub-cycles
        subCycle(GeometricField& gf, const label nSubCycles)
        :

            subCycleField<GeometricField>(gf),
            subCycleTime(const_cast<Time&>(gf.time()), nSubCycles)
        {}


    //- Destructor
    ~subCycle()
    {
        // End the subCycleTime, which restores the time state
        endSubCycle();
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
