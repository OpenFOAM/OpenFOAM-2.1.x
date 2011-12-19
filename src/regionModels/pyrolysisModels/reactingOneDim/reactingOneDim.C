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

\*---------------------------------------------------------------------------*/

#include "reactingOneDim.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactingOneDim, 0);

addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, mesh);
addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void reactingOneDim::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;
    time_.controlDict().lookup("maxDi") >> maxDiff_;

    coeffs().lookup("radFluxName") >> primaryRadFluxName_;
    coeffs().lookup("minimumDelta") >> minimumDelta_;
}


bool reactingOneDim::read()
{
    if (pyrolysisModel::read())
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


bool reactingOneDim::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


void reactingOneDim::updateQr()
{
    // Retrieve field from coupled region using mapped boundary conditions
    QrCoupled_.correctBoundaryConditions();

    // Update local Qr from coupled Qr field
    Qr_ == dimensionedScalar("zero", Qr_.dimensions(), 0.0);
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];

        scalarField& Qrp = Qr_.boundaryField()[patchI];

        // Qr is negative going out the solid
        // If the surface is emitting the radiative flux is set to zero
        Qrp = max(Qrp, scalar(0.0));
    }

    // Propagate Qr through 1-D regions
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];

        const scalarField& Qrp = Qr_.boundaryField()[patchI];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchI];

        forAll(Qrp, faceI)
        {
            const scalar Qr0 = Qrp[faceI];
            point Cf0 = Cf[faceI];
            const labelList& cells = boundaryFaceCells_[faceI];
            scalar kappaInt = 0.0;
            forAll(cells, k)
            {
                const label cellI = cells[k];
                const point& Cf1 = regionMesh().cellCentres()[cellI];
                const scalar delta = mag(Cf1 - Cf0);
                kappaInt += kappa_[cellI]*delta;
                Qr_[cellI] = Qr0*exp(-kappaInt);
                Cf0 = Cf1;
            }
        }
    }

    Qr_.correctBoundaryConditions();
}


void reactingOneDim::updatePhiGas()
{
    phiHsGas_ ==  dimensionedScalar("zero", phiHsGas_.dimensions(), 0.0);
    phiGas_ == dimensionedScalar("zero", phiGas_.dimensions(), 0.0);

    const speciesTable& gasTable = solidChemistry_->gasTable();

    forAll(gasTable, gasI)
    {
        tmp<volScalarField> tHsiGas = solidChemistry_->gasHs(T_, gasI);
        tmp<volScalarField> tRRiGas = solidChemistry_->RRg(gasI);

        const volScalarField& HsiGas = tHsiGas();
        const volScalarField& RRiGas = tRRiGas();

        const surfaceScalarField HsiGasf(fvc::interpolate(HsiGas));
        const surfaceScalarField RRiGasf(fvc::interpolate(RRiGas));

        forAll(intCoupledPatchIDs_, i)
        {
            const label patchI = intCoupledPatchIDs_[i];
            const scalarField& phiGasp = phiHsGas_.boundaryField()[patchI];

            forAll(phiGasp, faceI)
            {
                const labelList& cells = boundaryFaceCells_[faceI];
                scalar massInt = 0.0;
                forAllReverse(cells, k)
                {
                    const label cellI = cells[k];
                    massInt += RRiGas[cellI]*regionMesh().V()[cellI];
                    phiHsGas_[cellI] += massInt*HsiGas[cellI];
                }

                phiGas_.boundaryField()[patchI][faceI] += massInt;

                if (debug)
                {
                    Info<< " Gas : " << gasTable[gasI]
                        << " on patch : " << patchI
                        << " mass produced at face(local) : "
                        <<  faceI
                        << " is : " << massInt
                        << " [kg/s] " << endl;
                }
            }
        }
        tHsiGas().clear();
    }
}


void reactingOneDim::updateFields()
{
    updateQr();

    updatePhiGas();
}


void reactingOneDim::updateMesh(const scalarField& mass0)
{
    if (!moveMesh_)
    {
        return;
    }

    const scalarField newV(mass0/rho_);

    Info<< "Initial/final volumes = " << gSum(regionMesh().V()) << ", "
        << gSum(newV) << " [m3]" << endl;

    // move the mesh
    const labelList moveMap = moveMesh(regionMesh().V() - newV, minimumDelta_);

    // flag any cells that have not moved as non-reacting
    forAll(moveMap, i)
    {
        if (moveMap[i] == 0)
        {
            solidChemistry_->setCellReacting(i, false);
        }
    }
}


void reactingOneDim::solveContinuity()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveContinuity()" << endl;
    }

    solve
    (
        fvm::ddt(rho_)
     ==
      - solidChemistry_->RRg()
    );
}


void reactingOneDim::solveSpeciesMass()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveSpeciesMass()" << endl;
    }

    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho_, Yi)
         ==
            solidChemistry_->RRs(i)
        );

        if (moveMesh_)
        {
            surfaceScalarField phiRhoMesh
            (
                fvc::interpolate(Yi*rho_)*regionMesh().phi()
            );

            YiEqn -= fvc::div(phiRhoMesh);
        }

        YiEqn.solve(regionMesh().solver("Yi"));
        Yi.max(0.0);
        Yt += Yi;
    }

    Ys_[Ys_.size() - 1] = 1.0 - Yt;
}


void reactingOneDim::solveEnergy()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveEnergy()" << endl;
    }

    const volScalarField rhoCp(rho_*solidThermo_.Cp());

    const surfaceScalarField phiQr(fvc::interpolate(Qr_)*nMagSf());

    const surfaceScalarField phiGas(fvc::interpolate(phiHsGas_));

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T_)
      - fvm::laplacian(K_, T_)
     ==
        chemistrySh_
      + fvc::div(phiQr)
      + fvc::div(phiGas)
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

    Info<< "pyrolysis min/max(T) = " << min(T_).value() << ", "
        << max(T_).value() << endl;
}


void reactingOneDim::calculateMassTransfer()
{
    totalGasMassFlux_ = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        totalGasMassFlux_ += gSum(phiGas_.boundaryField()[patchI]);
    }

    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistrySh_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingOneDim::reactingOneDim(const word& modelType, const fvMesh& mesh)
:
    pyrolysisModel(modelType, mesh),
    solidChemistry_(solidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    kappa_(solidThermo_.kappa()),
    K_(solidThermo_.K()),
    rho_(solidThermo_.rho()),
    Ys_(solidThermo_.composition().Y()),
    T_(solidThermo_.T()),
    primaryRadFluxName_(coeffs().lookupOrDefault<word>("radFluxName", "Qr")),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    QrCoupled_
    (
        IOobject
        (
            primaryRadFluxName_,
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    Qr_
    (
        IOobject
        (
            "QrPyr",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        zeroGradientFvPatchVectorField::typeName
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0))
{
    if (active_)
    {
        read();
    }
}


reactingOneDim::reactingOneDim
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pyrolysisModel(modelType, mesh, dict),
    solidChemistry_(solidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    kappa_(solidThermo_.kappa()),
    K_(solidThermo_.K()),
    rho_(solidThermo_.rho()),
    Ys_(solidThermo_.composition().Y()),
    T_(solidThermo_.T()),
    primaryRadFluxName_(dict.lookupOrDefault<word>("radFluxName", "Qr")),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    QrCoupled_
    (
        IOobject
        (
            primaryRadFluxName_,
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    Qr_
    (
        IOobject
        (
            "QrPyr",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        zeroGradientFvPatchVectorField::typeName
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0))
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactingOneDim::~reactingOneDim()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar reactingOneDim::addMassSources(const label patchI, const label faceI)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchI)
        {
            index = i;
            break;
        }
    }

    const label localPatchId =  intCoupledPatchIDs_[index];

    const scalar massAdded = phiGas_.boundaryField()[localPatchId][faceI];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


scalar reactingOneDim::solidRegionDiffNo() const
{
    scalar DiNum = 0.0;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            regionMesh().surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(K_)
          / fvc::interpolate(Cp()*rho_)
        );

        DiNum = max(KrhoCpbyDelta.internalField())*time_.deltaTValue();
    }

    return DiNum;
}


scalar reactingOneDim::maxDiff() const
{
    return maxDiff_;
}


const volScalarField& reactingOneDim::rho() const
{
    return rho_;
}


const volScalarField& reactingOneDim::T() const
{
    return T_;
}


const tmp<volScalarField> reactingOneDim::Cp() const
{
    return solidThermo_.Cp();
}


const volScalarField& reactingOneDim::kappa() const
{
    return kappa_;
}


const volScalarField& reactingOneDim::K() const
{
    return K_;
}


const surfaceScalarField& reactingOneDim::phiGas() const
{
    return phiGas_;
}


void reactingOneDim::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();

    // Initialise all cells as able to react
    forAll(T_, cellI)
    {
        solidChemistry_->setCellReacting(cellI, true);
    }

    // De-activate reactions if pyrolysis region coupled to (valid) film
    if (filmCoupled_)
    {
        const volScalarField& filmDelta = filmDeltaPtr_();

        forAll(intCoupledPatchIDs_, i)
        {
            const label patchI = intCoupledPatchIDs_[i];
            const scalarField& filmDeltap = filmDelta.boundaryField()[patchI];

            forAll(filmDeltap, faceI)
            {
                const scalar filmDelta0 = filmDeltap[faceI];
                if (filmDelta0 > reactionDeltaMin_)
                {
                    const labelList& cells = boundaryFaceCells_[faceI];

                    // TODO: only limit cell adjacent to film?
                    //solidChemistry_->setCellNoReacting(cells[0])

                    // Propagate flag through 1-D region
                    forAll(cells, k)
                    {
                        solidChemistry_->setCellReacting(cells[k], false);
                    }
                }
            }
        }
    }
}


void reactingOneDim::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    const scalarField mass0 = rho_*regionMesh().V();

    solidChemistry_->solve
    (
        time().value() - time().deltaTValue(),
        time().deltaTValue()
    );

    solveContinuity();

    updateMesh(mass0);

    chemistrySh_ = solidChemistry_->Sh()();

    updateFields();

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }

    calculateMassTransfer();

    solidThermo_.correct();
}


void reactingOneDim::info() const
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << endl;

    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //
