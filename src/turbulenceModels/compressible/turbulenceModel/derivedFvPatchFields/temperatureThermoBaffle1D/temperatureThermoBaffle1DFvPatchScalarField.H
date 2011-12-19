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
    Foam::temperatureThermoBaffle1DFvPatchScalarField

Description
    Boundary which solves the 1D steady state heat transfer equation
    through a baffle.

SourceFiles
    temperatureThermoBaffle1DFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureThermoBaffle1DFvPatchScalarField_H
#define temperatureThermoBaffle1DFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
        Class temperatureThermoBaffle1DFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

template<class solidType>
class temperatureThermoBaffle1DFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

    class solidThermoData
    {
            // Solid thermo
            autoPtr<solidType> solidPtr_;

            // Solid dictionaries
            dictionary  transportDict_;
            dictionary  radiationDict_;
            dictionary  thermoDict_;
            dictionary  densityDict_;


    public:

        // Constructor from components
        solidThermoData(const dictionary& dict)
        :
            solidPtr_(new solidType(dict)),
            transportDict_(dict.subDict("transport")),
            radiationDict_(dict.subDict("radiation")),
            thermoDict_(dict.subDict("thermodynamics")),
            densityDict_(dict.subDict("density"))
        {}


        // Null constructor
        solidThermoData()
        :
            solidPtr_(),
            transportDict_(),
            radiationDict_(),
            thermoDict_(),
            densityDict_()
        {}


        // Destructor
        virtual ~solidThermoData()
        {}


        // Member Functions

            void write(Ostream& os) const
            {
                os.writeKeyword("transport");
                os << transportDict_  << nl;
                os.writeKeyword("radiation");
                os << radiationDict_ <<  nl;
                os.writeKeyword("thermodynamics");
                os << thermoDict_ << nl;
                os.writeKeyword("density");
                os << densityDict_ << nl;
            }


            scalar K(const scalar T) const
            {
                return solidPtr_().K(T);
            }
    };


        //- Name of the temperature field
        word TName_;

        //- Baffle is activated
        bool baffleActivated_;

        //- Baffle thickness [m]
        scalarField thickness_;

        //- Superficial heat source [W/m2]
        scalarField Qs_;

        //- Solid thermo
        autoPtr<solidThermoData> solid_;


public:

    //- Runtime type information
    TypeName("compressible::temperatureThermoBaffle1D");


    // Constructors

        //- Construct from patch and internal field
        temperatureThermoBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureThermoBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  temperatureThermoBaffle1DFvPatchScalarField onto a new patch
        temperatureThermoBaffle1DFvPatchScalarField
        (
            const temperatureThermoBaffle1DFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        temperatureThermoBaffle1DFvPatchScalarField
        (
            const temperatureThermoBaffle1DFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureThermoBaffle1DFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        temperatureThermoBaffle1DFvPatchScalarField
        (
            const temperatureThermoBaffle1DFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureThermoBaffle1DFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "temperatureThermoBaffle1DFvPatchScalarField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
