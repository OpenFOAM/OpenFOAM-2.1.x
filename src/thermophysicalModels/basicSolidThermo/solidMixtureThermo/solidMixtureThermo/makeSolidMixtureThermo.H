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

InClass
    Foam::solidMixtureThermo

\*---------------------------------------------------------------------------*/

#ifndef makeSolidMixtureThermo_H
#define makeSolidMixtureThermo_H

#include "addToRunTimeSelectionTable.H"
#include "solidMixtureThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolidMixtureThermo(CThermo,MixtureThermo,Mixture,Transport,Radiation,Thermo,Rho)\
                                                                              \
typedef MixtureThermo<Mixture<Transport<Radiation<Thermo<Rho> > > > >         \
    MixtureThermo##Mixture##Transport##Radiation##Thermo##Rho;                \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    MixtureThermo##Mixture##Transport##Radiation##Thermo##Rho,                \
    #MixtureThermo                                                            \
    "<"#Mixture"<"#Transport"<"#Radiation"<"#Thermo"<"#Rho">>>>>",            \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    CThermo,                                                                  \
    MixtureThermo##Mixture##Transport##Radiation##Thermo##Rho,                \
    mesh                                                                      \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    CThermo,                                                                  \
    MixtureThermo##Mixture##Transport##Radiation##Thermo##Rho,                \
    dictionary                                                                \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
