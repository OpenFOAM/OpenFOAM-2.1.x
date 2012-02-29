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
    Foam::solidMixtureThermo

Description
    Foam::solidMixtureThermo

SourceFiles
    solidMixtureThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidMixtureThermo_H
#define solidMixtureThermo_H

#include "basicSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solidMixtureThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class solidMixtureThermo
:
    public basicSolidThermo,
    public MixtureType
{
protected:

    // Protected data

        //- Thermal conductivity [W/m/K]
        volScalarField K_;


private:

    // Private Member Functions

        //- Calculate K
        void calculate();


        // Per patch helpers to get mixture properties on a patch

            tmp<scalarField> rho(const label patchI) const;

            tmp<scalarField> sigmaS(const label patchI) const;

            tmp<scalarField> kappa(const label patchI) const;

            tmp<scalarField> emissivity(const label patchI) const;

public:

    //- Runtime type information
    TypeName("solidMixtureThermo");


    // Constructors

        //- Construct from mesh
        solidMixtureThermo(const fvMesh&);

        //- Construct from mesh and dictionary
        solidMixtureThermo(const fvMesh&, const dictionary&);


    //- Destructor
    virtual ~solidMixtureThermo();


    // Member functions

        //- Update properties
        virtual void correct();

        //- Return the compostion of the solid mixture
        virtual MixtureType& composition()
        {
            return *this;
        }

        //- Return the compostion of the solid mixture
        virtual const MixtureType& composition() const
        {
            return *this;
        }


        // Derived thermal properties

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> K() const;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volSymmTensorField> directionalK() const
            {
                notImplemented("solidMixtureThermo::directionalK() const");
                return tmp<volSymmTensorField>(NULL);
            }

            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<volScalarField> Cp() const;

            //- Sensible enthalpy [J/(kg.K)]
            virtual tmp<volScalarField> hs() const;

            //- Heat of formation [J/kg]
            virtual tmp<volScalarField> Hf() const;


        // Per patch calculation

            //- Thermal conductivity [W/(m.K)]
            virtual tmp<scalarField> K(const label patchI) const;

            //- Thermal conductivity [W/(m.K)]
            virtual tmp<symmTensorField> directionalK(const label) const
            {
                notImplemented
                (
                    "solidMixtureThermo::directionalK(const label) const"
                );
                return tmp<symmTensorField>(NULL);
            }


            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<scalarField> Cp(const label patchI) const;

            //- Sensible enthalpy [J/(kg.K)]
            virtual tmp<scalarField> hs(const label patchI) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();

        //- Write the basicSolidThermo properties
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solidMixtureThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
