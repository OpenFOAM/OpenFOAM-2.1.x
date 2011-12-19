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
    Foam::multiComponentSolidMixture

Description
    Foam::multiComponentSolidMixture

\*---------------------------------------------------------------------------*/

#ifndef multiComponentSolidMixture_H
#define multiComponentSolidMixture_H

#include "PtrList.H"
#include "autoPtr.H"
#include "basicSolidMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class multiComponentSolidMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoSolidType>
class multiComponentSolidMixture
:
    public basicSolidMixture
{

    // Private data

        //- Solid data
        PtrList<ThermoSolidType> solidData_;


    // Private Member Functions

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Return molar fraction for component i in celli and at T
        scalar X(label i, label celli, scalar T) const;


public:


    // Constructors

        //- Construct from dictionary and mesh
        multiComponentSolidMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~multiComponentSolidMixture()
    {}


    // Member Functions

        //- Return the raw solid data
        const PtrList<ThermoSolidType>& solidData() const
        {
            return solidData_;
        }

        //- Read dictionary
        void read(const dictionary&);


        // Cell-based properties

            //- Density
            virtual scalar rho(scalar T, label celli) const;

            //- Absorption coefficient
            virtual scalar kappa(scalar T, label celli) const;

            //- Scatter coefficient
            virtual scalar sigmaS(scalar T, label celli) const;

            //- Thermal conductivity
            virtual scalar K(scalar T, label celli) const;

            //- Emissivity coefficient
            virtual scalar emissivity(scalar T, label celli) const;

            //- Formation enthalpy
            virtual scalar hf(scalar T, label celli) const;

            //- Sensible enthalpy
            virtual scalar hs(scalar T, label celli) const;

            //- Total enthalpy
            virtual scalar h(scalar T, label celli) const;

            //- Specific heat capacity
            virtual scalar Cp(scalar T, label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multiComponentSolidMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
