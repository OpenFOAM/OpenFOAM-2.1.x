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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::rhoChemistryModel&
Foam::combustionModels::rhoChemistryCombustionModel::pChemistry()
{
    return pChemistry_();
}

inline const Foam::rhoChemistryModel&
Foam::combustionModels::rhoChemistryCombustionModel::
pChemistry() const
{
    return pChemistry_();
}

inline Foam::tmp<Foam::volScalarField>
Foam::combustionModels::rhoChemistryCombustionModel::rho() const
{
    return pChemistry_->thermo().rho();
}

inline const Foam::hsReactionThermo&
Foam::combustionModels::rhoChemistryCombustionModel::thermo() const
{
    return pChemistry_->thermo();
}

inline Foam::hsReactionThermo&
Foam::combustionModels::rhoChemistryCombustionModel::thermo()
{
    return pChemistry_->thermo();
}

// ************************************************************************* //
