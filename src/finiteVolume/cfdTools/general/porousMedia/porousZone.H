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

Namespace
    Foam::porousMedia

Description
    Namespace for models related to porous media

Class
    Foam::porousZone

Description
    Porous zone definition based on cell zones.

    Porous zone definition based on cell zones and parameters obtained from a
    control dictionary constructed from the given stream. The orientation of
    the porous region is defined with the same notation as a coordinateSystem,
    but only a Cartesian coordinate system is valid.

    Implemented porosity models:

    powerLaw (\e C0 and \e C1 parameters)
    \f[
        S = - \rho C_0 |U|^{(C_1 - 1)} U
    \f]

    Darcy-Forchheimer (@e d and \e f parameters)
    \f[
        S = - (\mu \, d + \frac{\rho |U|}{2} \, f) U
    \f]


    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).

    The porousZones method porousZones::ddt() mirrors the normal fvm::ddt()
    method, but accounts for the effective volume of the cells.

    An example dictionary entry:
    \verbatim
    cat1
    {
        note            "some catalyst";
        coordinateSystem    system_10;
        porosity        0.809;
        intensity       0.001;      // optional
        mixingLength    0.0001;     // optional
        printCoeffs     yes;        // optional: feedback for the user

        Darcy
        {
            d   d [0 -2 0 0 0]  (-1000 -1000 5.3756e+07);
            f   f [0 -1 0 0 0]  (-1000 -1000 15.83);
        }
    }
    \endverbatim

See Also
    porousZones and coordinateSystems

SourceFiles
    porousZone.C
    porousZoneTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef porousZone_H
#define porousZone_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class porousZone Declaration
\*---------------------------------------------------------------------------*/

class porousZone
{
    // Private data

        //- Name of this zone, or a regular expression
        keyType key_;

        //- Reference to the finite volume mesh this zone is part of
        const fvMesh& mesh_;

        //- Dictionary containing the parameters
        dictionary dict_;

        //- Cell zone Ids
        labelList cellZoneIds_;

        //- Coordinate system used for the zone (Cartesian)
        coordinateSystem coordSys_;

        //- porosity of the zone (0 < porosity <= 1)
        //  Placeholder for treatment of temporal terms.
        //  Currently unused.
        scalar porosity_;

        //- Turbulent intensity as fraction of the velocity
        scalar intensity_;

        //- Turbulent length scale
        scalar mixingLength_;

        //- powerLaw coefficient C0
        scalar C0_;

        //- powerLaw coefficient C1
        scalar C1_;

        //- Darcy coefficient
        dimensionedTensor D_;

        //- Forchheimer coefficient
        dimensionedTensor F_;


    // Private Member Functions

        //- adjust negative resistance values to be multiplier of max value
        static void adjustNegativeResistance(dimensionedVector& resist);

        //- Power-law resistance
        template<class RhoFieldType>
        void addPowerLawResistance
        (
            scalarField& Udiag,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;


        //- Power-law resistance
        template<class RhoFieldType>
        void addPowerLawResistance
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;


        //- Disallow default bitwise copy construct
        porousZone(const porousZone&);

        //- Disallow default bitwise assignment
        void operator=(const porousZone&);


public:

    // Constructors

        //- Construct from components
        porousZone(const keyType& key, const fvMesh&, const dictionary&);

        //- Return clone
        autoPtr<porousZone> clone() const
        {
            notImplemented("autoPtr<porousZone> clone() const");
            return autoPtr<porousZone>(NULL);
        }

        //- Return pointer to new porousZone created on freestore from Istream
        class iNew
        {
            //- Reference to the finite volume mesh this zone is part of
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<porousZone> operator()(Istream& is) const
            {
                keyType key(is);
                dictionary dict(is);

                return autoPtr<porousZone>(new porousZone(key, mesh_, dict));
            }
        };


    //- Destructor
    virtual ~porousZone()
    {}


    // Member Functions

        // Access

            //- cellZone name
            const keyType& zoneName() const
            {
                return key_;
            }

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- cellZone numbers
            const labelList& zoneIds() const
            {
                return cellZoneIds_;
            }

            //- dictionary values used for the porousZone
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return coordinate system
            const coordinateSystem& coordSys() const
            {
                return coordSys_;
            }

            //- Return origin
            const point& origin() const
            {
                return coordSys_.origin();
            }

            //- Return axis
            vector axis() const
            {
                return coordSys_.axis();
            }

            //- Return porosity
            scalar porosity() const
            {
                return porosity_;
            }

            //- Edit access to porosity
            scalar& porosity()
            {
                return porosity_;
            }

            //- Return turbulent intensity
            scalar intensity() const
            {
                return intensity_;
            }

            //- Edit access to turbulent intensity
            scalar& intensity()
            {
                return intensity_;
            }

            //- Return turbulent length scale
            scalar mixingLength() const
            {
                return mixingLength_;
            }

            //- Edit access to turbulent length scale
            scalar& mixingLength()
            {
                return mixingLength_;
            }


        //- Modify time derivative elements according to porosity
        template<class Type>
        void modifyDdt(fvMatrix<Type>&) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        void addResistance(fvVectorMatrix& UEqn) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation giving rho and mu.
        void addResistance
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal.
        //  Optionally correct the processor BCs of AU.
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC = true
        ) const;

        //- Write the porousZone dictionary
        virtual void writeDict(Ostream&, bool subDict = true) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const porousZone&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porousZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
