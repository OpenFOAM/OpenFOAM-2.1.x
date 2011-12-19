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

#include "GidaspowErgunWenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(GidaspowErgunWenYu, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        GidaspowErgunWenYu,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::GidaspowErgunWenYu::GidaspowErgunWenYu
(
    const dictionary& interfaceDict,
    const volScalarField& alpha,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    dragModel(interfaceDict, alpha, phase1, phase2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::GidaspowErgunWenYu::~GidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::GidaspowErgunWenYu::K
(
    const volScalarField& Ur
) const
{
    volScalarField beta(max(scalar(1) - alpha_, scalar(1.0e-6)));
    volScalarField d(phase1_.d());
    volScalarField bp(pow(beta, -2.65));
    volScalarField Re(max(Ur*d/phase2_.nu(), scalar(1.0e-3)));

    volScalarField Cds(24.0*(1.0 + 0.15*pow(Re, 0.687))/Re);

    forAll(Re, celli)
    {
        if (Re[celli] > 1000.0)
        {
            Cds[celli] = 0.44;
        }
    }

    // Wen and Yu (1966)
    tmp<volScalarField> tKWenYu(0.75*Cds*phase2_.rho()*Ur*bp/d);
    volScalarField& KWenYu = tKWenYu();

    // Ergun
    forAll (beta, cellj)
    {
        if (beta[cellj] <= 0.8)
        {
            KWenYu[cellj] =
                150.0*alpha_[cellj]*phase2_.nu().value()*phase2_.rho().value()
               /sqr(beta[cellj]*d[cellj])
              + 1.75*phase2_.rho().value()*Ur[cellj]
               /(beta[cellj]*d[cellj]);
        }
    }

    return tKWenYu;
}


// ************************************************************************* //
