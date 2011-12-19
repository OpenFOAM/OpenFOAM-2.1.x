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
    Foam::thermalPorousZone

Description
    Porous zone definition based on cell zones including terms for energy
    equations.

See Also
    porousZone, thermalPorousZones and coordinateSystems with run-time
    selectable thermal model

SourceFiles
    thermalPorousZone.C

\*---------------------------------------------------------------------------*/

#ifndef thermalPorousZone_H
#define thermalPorousZone_H

#include "porousZone.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class basicThermo;

/*---------------------------------------------------------------------------*\
                       Class thermalPorousZone Declaration
\*---------------------------------------------------------------------------*/

class thermalPorousZone
:
    public porousZone
{
    // Private data

        //- Disallow default bitwise copy construct
        thermalPorousZone(const thermalPorousZone&);

        //- Disallow default bitwise assignment
        void operator=(const thermalPorousZone&);

        //- Thermal model
        autoPtr<porousMedia::thermalModel> model_;


public:

    // Constructors

        //- Construct from components
        thermalPorousZone
        (
            const keyType& key,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Return clone
        autoPtr<thermalPorousZone> clone() const
        {
            notImplemented("autoPtr<thermalPorousZone> clone() const");
            return autoPtr<thermalPorousZone>(NULL);
        }

        //- Return pointer to new thermalPorousZone
        //  created on freestore from Istream
        class iNew
        {
            //- Reference to the finite volume mesh this zone is part of
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<thermalPorousZone> operator()(Istream& is) const
            {
                keyType key(is);
                dictionary dict(is);

                return autoPtr<thermalPorousZone>
                (
                    new thermalPorousZone(key, mesh_, dict)
                );
            }
        };


    //- Destructor
    virtual ~thermalPorousZone()
    {}


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
