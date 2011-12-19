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

#include "sigWriteNow.H"
#include "error.H"
#include "JobInfo.H"
#include "IOstreams.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Signal number to catch
int Foam::sigWriteNow::signal_
(
    debug::optimisationSwitch("writeNowSignal", -1)
);

static Foam::Time* runTimePtr_ = NULL;


struct sigaction Foam::sigWriteNow::oldAction_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigWriteNow::sigHandler(int)
{
    Info<< "sigWriteNow :"
        << " setting up write at end of the next iteration" << nl << endl;
    runTimePtr_->writeOnce();

    //// Throw signal (to old handler)
    //raise(signal_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigWriteNow::sigWriteNow()
{}


Foam::sigWriteNow::sigWriteNow(const bool verbose, Time& runTime)
{
    if (signal_ >= 0)
    {
        // Store runTime
        runTimePtr_ = &runTime;

        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(signal_, &newAction, &oldAction_) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigWriteNow::sigWriteNow(const bool, const Time&)"
            )   << "Cannot set " << signal_ << " trapping"
                << abort(FatalError);
        }

        if (verbose)
        {
            Info<< "sigWriteNow :"
                << " Enabling writing upon signal " << signal_
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigWriteNow::~sigWriteNow()
{
    // Reset old handling
    if (signal_ > 0)
    {
        if (sigaction(signal_, &oldAction_, NULL) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigWriteNow::~sigWriteNow()"
            )   << "Cannot reset " << signal_ << " trapping"
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sigWriteNow::active() const
{
    return signal_ > 0;
}


// ************************************************************************* //
