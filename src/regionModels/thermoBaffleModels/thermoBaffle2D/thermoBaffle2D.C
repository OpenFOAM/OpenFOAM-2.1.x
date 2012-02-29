/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "thermoBaffle2D.H"

#include "fvm.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermoBaffleModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoBaffle2D, 0);

addToRunTimeSelectionTable(thermoBaffleModel, thermoBaffle2D, mesh);
addToRunTimeSelectionTable(thermoBaffleModel, thermoBaffle2D, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


bool thermoBaffle2D::read()
{
    this->solution().lookup("nNonOrthCorr") >> nNonOrthCorr_;
    return regionModel1D::read();
}


bool thermoBaffle2D::read(const dictionary& dict)
{
    this->solution().lookup("nNonOrthCorr") >> nNonOrthCorr_;
    return regionModel1D::read(dict);
}


void thermoBaffle2D::solveEnergy()
{
    if (debug)
    {
        Info<< "thermoBaffle2D::solveEnergy()" << endl;
    }

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    tmp<volScalarField> tQ
    (
        new volScalarField
        (
            IOobject
            (
                "tQ",
                regionMesh().time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    volScalarField& Q = tQ();

    volScalarField rhoCp("rhoCp", thermo_->rho()*thermo_->Cp()());
    volScalarField K("K", thermo_->K());


    //If region is one-dimension variable thickness
    if (oneD_ && !constantThickness_)
    {
        // Scale K and rhoCp and fill Q in the internal baffle region.
        const label patchI = intCoupledPatchIDs_[0];
        const polyPatch& ppCoupled = rbm[patchI];

        forAll(ppCoupled, localFaceI)
        {
            const labelList& cells = boundaryFaceCells_[localFaceI];
            forAll (cells, i)
            {
                const label cellId = cells[i];

                Q[cellId] =
                    Qs_.boundaryField()[patchI][localFaceI]
                    /thickness_[localFaceI];

                rhoCp[cellId] *= delta_.value()/thickness_[localFaceI];

                K[cellId] *= delta_.value()/thickness_[localFaceI];
            }
        }
    }
    else
    {
        Q = Q_;
    }

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T_)
      - fvm::laplacian(K, T_)
     ==
        Q
    );

    if (moveMesh_)
    {
        surfaceScalarField phiMesh
        (
            fvc::interpolate(rhoCp*T_)*regionMesh().phi()
        );

        TEqn -= fvc::div(phiMesh);
    }

    TEqn.relax();
    TEqn.solve();

    Info<< "T gas min/max   = " << min(T_).value() << ", "
        << max(T_).value() << endl;

    thermo_->correct();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


thermoBaffle2D::thermoBaffle2D
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    thermoBaffleModel(modelType, mesh, dict),
    nNonOrthCorr_(readLabel(solution().lookup("nNonOrthCorr"))),
    thermo_(basicSolidThermo::New(regionMesh(), dict)),
    T_(thermo_->T()),
    Qs_
    (
        IOobject
        (
            "Qs",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimArea/dimTime,
            pTraits<scalar>::zero
        )
    ),
    Q_
    (
        IOobject
        (
            "Q",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimVolume/dimTime,
            pTraits<scalar>::zero
        )
    )
{
    init();
    thermo_->correct();
}


thermoBaffle2D::thermoBaffle2D
(
    const word& modelType,
    const fvMesh& mesh
)
:
    thermoBaffleModel(modelType, mesh),
    nNonOrthCorr_(readLabel(solution().lookup("nNonOrthCorr"))),
    thermo_(basicSolidThermo::New(regionMesh())),
    T_(thermo_->T()),
    Qs_
    (
        IOobject
        (
            "Qs",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimArea/dimTime,
            pTraits<scalar>::zero
        )
    ),
    Q_
    (
        IOobject
        (
            "Q",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar
        (
            "zero",
            dimEnergy/dimVolume/dimTime,
            pTraits<scalar>::zero
        )
    )
{
    init();
    thermo_->correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoBaffle2D::~thermoBaffle2D()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoBaffle2D::init()
{
    if (oneD_ && !constantThickness_)
    {
        label patchI = intCoupledPatchIDs_[0];
        const label Qsb = Qs_.boundaryField()[patchI].size();
        if (Qsb!= thickness_.size())
        {
            FatalErrorIn
            (
                "thermoBaffle2D::thermoBaffle2D"
                "("
                "   const word& modelType,"
                "   const fvMesh& mesh,"
                "   const dictionary& dict"
                ")"
            )   << "the boundary field of Qs is "
                << Qsb << " and " << nl
                << "the field 'thickness' is " << thickness_.size() << nl
                << exit(FatalError);
        }
    }
}


void thermoBaffle2D::preEvolveRegion()
{}


void thermoBaffle2D::evolveRegion()
{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }
}


const tmp<volScalarField> thermoBaffle2D::Cp() const
{
    return thermo_->Cp();
}


const volScalarField& thermoBaffle2D::kappa() const
{
    return thermo_->kappa();
}


const volScalarField& thermoBaffle2D::rho() const
{
    return thermo_->rho();
}


const volScalarField& thermoBaffle2D::K() const
{
    return thermo_->K();
}


const volScalarField& thermoBaffle2D::T() const
{
    return T_;
}


const basicSolidThermo& thermoBaffle2D::thermo() const
{
    return thermo_;
}


void thermoBaffle2D::info() const
{
    Info<< indent << "min/max(T) = " << min(T_).value() << ", "
        << max(T_).value() << nl;

    const labelList& coupledPatches = intCoupledPatchIDs();
    forAll (coupledPatches, i)
    {
        const label patchI = coupledPatches[i];
        const fvPatchScalarField& pT = T_.boundaryField()[patchI];
        const word patchName = regionMesh().boundary()[patchI].name();
        Info << indent << "Q : " << patchName << indent <<
            gSum
            (
                mag(regionMesh().Sf().boundaryField()[patchI])
              * pT.snGrad()
              * thermo_->K(patchI)
            ) << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace thermoBaffleModels
} // end namespace regionModels
} // end namespace Foam

// ************************************************************************* //
