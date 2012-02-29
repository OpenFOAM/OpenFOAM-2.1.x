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

\*---------------------------------------------------------------------------*/

#include "solidMixtureThermo.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
void Foam::solidMixtureThermo<MixtureType>::calculate()
{

    scalarField& rhoCells = rho_.internalField();
    scalarField& KCells = K_.internalField();
    scalarField& kappaCells = kappa_.internalField();
    scalarField& sigmaSCells = sigmaS_.internalField();
    scalarField& emissivityCells = emissivity_.internalField();

    const volScalarField::InternalField& iT = T_.internalField();

    forAll(iT, celli)
    {
        rhoCells[celli] = MixtureType::rho(iT[celli], celli);
        kappaCells[celli] = MixtureType::kappa(iT[celli], celli);
        sigmaSCells[celli] = MixtureType::sigmaS(iT[celli], celli);
        KCells[celli] = MixtureType::K(iT[celli], celli);
        emissivityCells[celli] = MixtureType::emissivity(iT[celli], celli);
    }

    const volScalarField::GeometricBoundaryField& bT = T_.boundaryField();

    forAll(bT, patchI)
    {
        rho_.boundaryField()[patchI] == this->rho(patchI)();
        K_.boundaryField()[patchI] == this->K(patchI)();
        kappa_.boundaryField()[patchI] == this->kappa(patchI)();
        sigmaS_.boundaryField()[patchI] == this->sigmaS(patchI)();
        emissivity_.boundaryField()[patchI] == this->emissivity(patchI)();
    }
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::rho
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
    const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> tRho(new scalarField(patchT.size()));
    scalarField& Rho = tRho();

    forAll(patchT, celli)
    {
        Rho[celli] = MixtureType::rho(patchT[celli], cells[celli]);
    }

    return tRho;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::sigmaS
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
    const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> tsigmaS(new scalarField(patchT.size()));
    scalarField& sigmaS = tsigmaS();

    forAll(patchT, celli)
    {
        sigmaS[celli] =
            MixtureType::sigmaS(patchT[celli], cells[celli]);
    }

    return tsigmaS;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::kappa
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
   const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> tKappa(new scalarField(patchT.size()));
    scalarField& kappa = tKappa();

    forAll(patchT, celli)
    {
        kappa[celli] =
            MixtureType::kappa(patchT[celli], cells[celli]);
    }

    return tKappa;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::emissivity
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
    const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> te(new scalarField(patchT.size()));
    scalarField& e = te();

    forAll(patchT, celli)
    {
        e[celli] = MixtureType::emissivity(patchT[celli], cells[celli]);
    }

    return te;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::solidMixtureThermo<MixtureType>::solidMixtureThermo
(
    const fvMesh& mesh
)
:
    basicSolidThermo(mesh),
    MixtureType(*this, mesh),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    )
{
    calculate();
}


template<class MixtureType>
Foam::solidMixtureThermo<MixtureType>::solidMixtureThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicSolidThermo(mesh, dict),
    MixtureType(*this, mesh),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    )
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::solidMixtureThermo<MixtureType>::~solidMixtureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::solidMixtureThermo<MixtureType>::correct()
{
    calculate();
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::solidMixtureThermo<MixtureType>::K() const
{
    return K_;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::solidMixtureThermo<MixtureType>::Cp() const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/(dimMass*dimTemperature)
        )
    );
    volScalarField& Cp = tCp();

    const volScalarField::InternalField& iT = T_.internalField();

    forAll(iT, celli)
    {
        Cp[celli] = MixtureType::Cp(T_[celli], celli);
    }

    volScalarField::GeometricBoundaryField& bCp = Cp.boundaryField();

    forAll(bCp, patchI)
    {
        bCp[patchI] == this->Cp(patchI)();
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::solidMixtureThermo<MixtureType>::hs() const
{
    tmp<volScalarField> ths
    (
        new volScalarField
        (
            IOobject
            (
                "Hs",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/(dimMass*dimTemperature)
        )
    );
    volScalarField& hs = ths();

    const volScalarField::InternalField& iT = T_.internalField();

    forAll(iT, celli)
    {
        hs[celli] = MixtureType::hs(T_[celli], celli);
    }

    volScalarField::GeometricBoundaryField& bHs = hs.boundaryField();

    forAll(bHs, patchI)
    {
        bHs[patchI] == this->hs(patchI)();
    }

    return ths;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::solidMixtureThermo<MixtureType>::Hf() const
{
    tmp<volScalarField> thF
    (
        new volScalarField
        (
            IOobject
            (
                "hF",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/(dimMass*dimTemperature)
        )
    );
    volScalarField& hf = thF();

    const volScalarField::InternalField& iT = T_.internalField();

    forAll(iT, celli)
    {
        hf[celli] = MixtureType::hf(T_[celli], celli);
    }

    volScalarField::GeometricBoundaryField& bhf = hf.boundaryField();

    forAll(bhf, patchI)
    {
        bhf[patchI] == this->Hf(patchI)();
    }

    return thF;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::K
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
    const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> tK(new scalarField(patchT.size()));
    scalarField& K = tK();

    forAll(patchT, celli)
    {
        K[celli] = MixtureType::K(patchT[celli], cells[celli]);
    }

    return tK;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::Cp
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
    const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> tCp(new scalarField(patchT.size()));
    scalarField& Cp = tCp();

    forAll(patchT, celli)
    {
        Cp[celli] = MixtureType::Cp(patchT[celli], cells[celli]);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::hs
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
    const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> ths(new scalarField(patchT.size()));
    scalarField& hs = ths();

    forAll(patchT, celli)
    {
        hs[celli] = MixtureType::hs(patchT[celli], cells[celli]);
    }

    return ths;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::solidMixtureThermo<MixtureType>::Hf
(
    const label patchI
) const
{
    const scalarField& patchT = T_.boundaryField()[patchI];
    const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    const labelUList& cells = pp.faceCells();

    tmp<scalarField> tHf(new scalarField(patchT.size()));
    scalarField& Hf = tHf();

    forAll(patchT, celli)
    {
        Hf[celli] = MixtureType::hf(patchT[celli], cells[celli]);
    }

    return tHf;
}


template<class MixtureType>
bool Foam::solidMixtureThermo<MixtureType>::read()
{
    if (basicSolidThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


template<class MixtureType>
bool Foam::solidMixtureThermo<MixtureType>::writeData(Ostream& os) const
{
     bool ok = basicSolidThermo::writeData(os);
     return ok && os.good();
}


// ************************************************************************* //
