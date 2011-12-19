/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

#include "singleStepCombustion.H"
#include "fvmSup.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
singleStepCombustion<CombThermoType, ThermoType>
::singleStepCombustion
(
    const word& modelType, const fvMesh& mesh
)
:
    CombThermoType(modelType, mesh),
    singleMixture_
    (
        dynamic_cast<singleStepReactingMixture<ThermoType>&>(this->thermo())
    ),
    wFuel_
    (
        IOobject
        (
            "wFuel",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
singleStepCombustion<CombThermoType, ThermoType>
::~singleStepCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::fvScalarMatrix>
singleStepCombustion<CombThermoType, ThermoType>::R
(
    const volScalarField& Y
) const
{
    const label specieI = this->thermo_->composition().species()[Y.name()];

    const volScalarField wSpecie
    (
        wFuel_*singleMixture_.specieStoichCoeffs()[specieI]
    );

    return wSpecie + fvm::Sp(0.0*wSpecie, Y);
}


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
singleStepCombustion< CombThermoType, ThermoType>::Sh() const
{
    const label fuelI = singleMixture_.fuelIndex();
    const volScalarField& YFuel = this->thermo_->composition().Y(fuelI);

    return -singleMixture_.qFuel()*(R(YFuel) & YFuel);
}


template<class CombThermoType, class ThermoType>
Foam::tmp<Foam::volScalarField>
singleStepCombustion< CombThermoType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh().V()*Sh()();
    }
    return tdQ;
}


template<class CombThermoType, class ThermoType>
bool singleStepCombustion< CombThermoType, ThermoType>::read()
{
    if (CombThermoType::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
