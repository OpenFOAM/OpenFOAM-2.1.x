/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "ThermoCloud.H"
#include "ThermoParcel.H"

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setModels()
{
    heatTransferModel_.reset
    (
        HeatTransferModel<ThermoCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    TIntegrator_.reset
    (
        scalarIntegrationScheme::New
        (
            "T",
            this->solution().integrationSchemes()
        ).ptr()
    );

    this->subModelProperties().lookup("radiation") >> radiation_;
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::cloudReset(ThermoCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    heatTransferModel_.reset(c.heatTransferModel_.ptr());
    TIntegrator_.reset(c.TIntegrator_.ptr());

    radiation_ = c.radiation_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        rho,
        U,
        thermo.thermo().mu(),
        g,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(NULL),
    constProps_(this->particleProperties(), this->solution().active()),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p()),
    heatTransferModel_(NULL),
    TIntegrator_(NULL),
    radiation_(false),
    hsTrans_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy, 0.0)
        )
    ),
    hsCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTemperature, 0.0)
        )
    )
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    ThermoCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermoCloud(),
    cloudCopyPtr_(NULL),
    constProps_(c.constProps_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(c.heatTransferModel_->clone()),
    TIntegrator_(c.TIntegrator_->clone()),
    radiation_(c.radiation_),
    hsTrans_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsTrans()
        )
    ),
    hsCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "hsCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            c.hsCoeff()
        )
    )
{}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const ThermoCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(NULL),
    constProps_(),
    thermo_(c.thermo()),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(NULL),
    TIntegrator_(NULL),
    radiation_(false),
    hsTrans_(NULL),
    hsCoeff_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoCloud<CloudType>::~ThermoCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.T() = constProps_.T0();
    parcel.Cp() = constProps_.Cp0();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ThermoCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    hsTrans_->field() = 0.0;
    hsCoeff_->field() = 0.0;
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::relaxSources
(
    const ThermoCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    this->relax(hsTrans_(), cloudOldTime.hsTrans(), "hs");
    this->relax(hsCoeff_(), cloudOldTime.hsCoeff(), "hs");
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(hsTrans_(), "hs");
    this->scale(hsCoeff_(), "hs");
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::preEvolve()
{
    CloudType::preEvolve();

    this->pAmbient() = thermo_.thermo().p().average().value();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<ThermoCloud<CloudType> > td(*this);

        this->solve(td);
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    typedef typename particle::TrackingData<ThermoCloud<CloudType> > tdType;

    tdType td(*this);

    Cloud<parcelType>::template autoMap<tdType>(td, mapper);
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::info()
{
    CloudType::info();

    Info<< "    Temperature min/max             = " << Tmin() << ", " << Tmax()
        << endl;
}


// ************************************************************************* //
