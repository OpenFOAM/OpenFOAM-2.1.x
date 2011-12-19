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
    Foam::egrMixture

Description
    Foam::egrMixture

SourceFiles
    egrMixture.C

\*---------------------------------------------------------------------------*/

#ifndef egrMixture_H
#define egrMixture_H

#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class egrMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class egrMixture
:
    public basicMultiComponentMixture
{
    // Private data

        static const int nSpecies_ = 3;
        static const char* specieNames_[3];

        dimensionedScalar stoicRatio_;

        ThermoType fuel_;
        ThermoType oxidant_;
        ThermoType products_;

        mutable ThermoType mixture_;

        //- Mixture fraction
        volScalarField& ft_;

        //- Regress variable
        volScalarField& b_;

        //- Residual gases
        volScalarField& egr_;

        //- Construct as copy (not implemented)
        egrMixture(const egrMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        egrMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~egrMixture()
    {}


    // Member functions

        const dimensionedScalar& stoicRatio() const
        {
            return stoicRatio_;
        }

        const ThermoType& mixture
        (
            const scalar,
            const scalar,
            const scalar
        ) const;

        const ThermoType& cellMixture(const label celli) const
        {
            return mixture(ft_[celli], b_[celli], egr_[celli]);
        }

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                b_.boundaryField()[patchi][facei],
                egr_.boundaryField()[patchi][facei]
            );
        }

        const ThermoType& cellReactants(const label celli) const
        {
            return mixture(ft_[celli], 1, egr_[celli]);
        }

        const ThermoType& patchFaceReactants
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                1,
                egr_.boundaryField()[patchi][facei]
            );
        }

        const ThermoType& cellProducts(const label celli) const
        {
            return mixture(ft_[celli], 0, 0);
        }

        const ThermoType& patchFaceProducts
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                0,
                0
            );
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        const ThermoType& getLocalThermo(const label specieI) const;


        // Per specie properties

            //- Number of moles []
            virtual scalar nMoles(const label specieI) const;

            //- Molecular weight [kg/kmol]
            virtual scalar W(const label specieI) const;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp(const label specieI, const scalar T) const;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv(const label specieI, const scalar T) const;

            //- Enthalpy [J/kg]
            virtual scalar H(const label specieI, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs(const label specieI, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label specieI) const;

            //- Entropy [J/(kg K)]
            virtual scalar S(const label specieI, const scalar T) const;

            //- Internal energy [J/kg]
            virtual scalar E(const label specieI, const scalar T) const;

            //- Gibbs free energy [J/kg]
            virtual scalar G(const label specieI, const scalar T) const;

            //- Helmholtz free energy [J/kg]
            virtual scalar A(const label specieI, const scalar T) const;


        // Per specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu(const label specieI, const scalar T) const;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa(const label specieI, const scalar T) const;

            //- Thermal diffusivity [kg/m/s]
            virtual scalar alpha(const label specieI, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "egrMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
