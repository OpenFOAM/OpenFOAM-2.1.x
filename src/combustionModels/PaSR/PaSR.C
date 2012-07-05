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

#include "PaSR.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType>
Foam::combustionModels::PaSR<CombThermoType>::PaSR
(
    const word& modelType,
    const fvMesh& mesh
)
:
    CombThermoType(modelType, mesh),
    Cmix_(this->coeffs().lookup("Cmix")),
    turbulentReaction_(this->coeffs().lookup("turbulentReaction")),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("kappa", dimless, 0.0)
    ),
    useReactionRate_(this->coeffs().lookupOrDefault("useReactionRate", false))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType>
Foam::combustionModels::PaSR<CombThermoType>::~PaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<class CombThermoType>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<CombThermoType>::tc() const
{
    return this->pChemistry_->tc();
}


template<class CombThermoType>
void Foam::combustionModels::PaSR<CombThermoType>::correct()
{
    if (this->active())
    {
        if (!useReactionRate_)
        {
            this->pChemistry_->solve
            (
                this->mesh().time().value()-this->mesh().time().deltaTValue(),
                this->mesh().time().deltaTValue()
            );
        }
        else
        {
            this->pChemistry_->calculate();
        }

        if (turbulentReaction_)
        {
            tmp<volScalarField> trho(this->rho());
            const volScalarField& rho = trho();
            tmp<volScalarField> tepsilon(this->turbulence().epsilon());
            const volScalarField& epsilon = tepsilon();
            tmp<volScalarField> tmuEff(this->turbulence().muEff());
            const volScalarField& muEff = tmuEff();
            tmp<volScalarField> ttc(tc());
            const volScalarField& tc = ttc();

            const scalar dt = this->mesh().time().deltaTValue();

            forAll(epsilon, i)
            {
                if (epsilon[i] > 0)
                {
                    scalar tk =
                        Cmix_.value()
                       *Foam::sqrt(muEff[i]/rho[i]/(epsilon[i] + SMALL));

                    // Chalmers PaSR model
                    if (!useReactionRate_)
                    {
                        kappa_[i] = (dt + tc[i])/(dt + tc[i] + tk);
                    }
                    else
                    {
                        kappa_[i] = tc[i]/(tc[i] + tk);
                    }
                }
                else
                {
                    // Return to laminar combustion
                    kappa_[i] = 1.0;
                }
            }
        }
        else
        {
            kappa_ = 1.0;
        }
    }
}


template<class CombThermoType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR<CombThermoType>::R(const volScalarField& Y) const
{

    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimMass/dimTime)
    );

    fvScalarMatrix& Su = tSu();

    if (this->active())
    {
        const label specieI = this->thermo().composition().species()[Y.name()];

        Su += kappa_*this->pChemistry_->RR(specieI);
    }

    return tSu;
}


template<class CombThermoType>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<CombThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        volScalarField& dQ = tdQ();
        dQ = kappa_*this->pChemistry_->dQ();
    }

    return tdQ;
}


template<class CombThermoType>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR<CombThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        scalarField& Sh = tSh();
        Sh = kappa_*this->pChemistry_->Sh();
    }

    return tSh;
}


template<class CombThermoType>
bool Foam::combustionModels::PaSR<CombThermoType>::read()
{
    if (CombThermoType::read())
    {
        this->coeffs().lookup("Cmix") >> Cmix_;
        this->coeffs().lookup("turbulentReaction") >> turbulentReaction_;
        this->coeffs().lookup("useReactionRate") >> useReactionRate_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
