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

#include "multiComponentSolidMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoSolidType>
void Foam::multiComponentSolidMixture<ThermoSolidType>::correctMassFractions()
{
    volScalarField Yt("Yt", Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }


}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::X
(
    label iComp, label celli, scalar T
) const
{
    scalar rhoInv = 0.0;
    forAll(solidData_, i)
    {
        rhoInv += Y_[i][celli]/solidData_[i].rho(T);
    }

    scalar X = Y_[iComp][celli]/solidData_[iComp].rho(T);

    return (X/rhoInv);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoSolidType>
Foam::multiComponentSolidMixture<ThermoSolidType>::multiComponentSolidMixture
(
    const dictionary& thermoSolidDict,
    const fvMesh& mesh
)
:
    basicSolidMixture
    (
        thermoSolidDict.lookup("solidComponents"),
        mesh
    ),
    solidData_(components_.size())
{

    forAll(components_, i)
    {
        solidData_.set
        (
            i,
            new ThermoSolidType
            (
                thermoSolidDict.subDict(components_[i] + "Coeffs")
            )
        );
    }
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::rho
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].rho(T)*X(i, celli, T);
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::hf
(
    scalar, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].hf()*Y_[i][celli];
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::hs
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].hs(T)*Y_[i][celli];
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::h
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].h(T)*Y_[i][celli];
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::kappa
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].kappa(T)*X(i, celli, T);
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::sigmaS
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].sigmaS(T)*X(i, celli, T);
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::K
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].K(T)*X(i, celli, T);
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::emissivity
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].emissivity(T)*X(i, celli, T);
    }
    return tmp;
}


template<class ThermoSolidType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoSolidType>::Cp
(
    scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].Cp(T)*Y_[i][celli];
    }
    return tmp;
}


template<class ThermoSolidType>
void Foam::multiComponentSolidMixture<ThermoSolidType>::read
(
    const dictionary& thermoDict
)
{
    forAll(components_, i)
    {
        solidData_[i] =
            ThermoSolidType(thermoDict.subDict(components_[i] + "Coeffs"));
    }
}


// ************************************************************************* //
