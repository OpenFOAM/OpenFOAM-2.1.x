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

#include "mpi.h"

#include "UPstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "PstreamGlobals.H"
#include "SubList.H"

#include <cstring>
#include <cstdlib>
#include <csignal>

#if defined(WM_SP)
#   define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
#   define MPI_SCALAR MPI_DOUBLE
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("GAMMANP", "number of instances");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::init(int& argc, char**& argv)
{
    MPI_Init(&argc, &argv);

    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcNo_);

    if (debug)
    {
        Pout<< "UPstream::init : initialised with numProcs:" << numprocs
            << " myProcNo:" << myProcNo_ << endl;
    }

    if (numprocs <= 1)
    {
        FatalErrorIn("UPstream::init(int& argc, char**& argv)")
            << "bool IPstream::init(int& argc, char**& argv) : "
               "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    procIDs_.setSize(numprocs);

    forAll(procIDs_, procNo)
    {
        procIDs_[procNo] = procNo;
    }

    setParRun();

#   ifndef SGIMPI
    string bufferSizeName = getEnv("MPI_BUFFER_SIZE");

    if (bufferSizeName.size())
    {
        int bufferSize = atoi(bufferSizeName.c_str());

        if (bufferSize)
        {
            MPI_Buffer_attach(new char[bufferSize], bufferSize);
        }
    }
    else
    {
        FatalErrorIn("UPstream::init(int& argc, char**& argv)")
            << "UPstream::init(int& argc, char**& argv) : "
            << "environment variable MPI_BUFFER_SIZE not defined"
            << Foam::abort(FatalError);
    }
#   endif

    int processorNameLen;
    char processorName[MPI_MAX_PROCESSOR_NAME];

    MPI_Get_processor_name(processorName, &processorNameLen);

    //signal(SIGABRT, stop);

    // Now that nprocs is known construct communication tables.
    initCommunicationSchedule();

    return true;
}


void Foam::UPstream::exit(int errnum)
{
    if (debug)
    {
        Pout<< "UPstream::exit." << endl;
    }

#   ifndef SGIMPI
    int size;
    char* buff;
    MPI_Buffer_detach(&buff, &size);
    delete[] buff;
#   endif

    if (PstreamGlobals::outstandingRequests_.size())
    {
        label n = PstreamGlobals::outstandingRequests_.size();
        PstreamGlobals::outstandingRequests_.clear();

        WarningIn("UPstream::exit(int)")
            << "There are still " << n << " outstanding MPI_Requests." << endl
            << "This means that your code exited before doing a"
            << " UPstream::waitRequests()." << endl
            << "This should not happen for a normal code exit."
            << endl;
    }

    if (errnum == 0)
    {
        MPI_Finalize();
        ::exit(errnum);
    }
    else
    {
        MPI_Abort(MPI_COMM_WORLD, errnum);
    }
}


void Foam::UPstream::abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}


void Foam::reduce(scalar& Value, const sumOp<scalar>& bop, const int tag)
{
    if (Pstream::debug)
    {
        Pout<< "Foam::reduce : value:" << Value << endl;
    }

    if (!UPstream::parRun())
    {
        return;
    }

    if (UPstream::nProcs() <= UPstream::nProcsSimpleSum)
    {
        if (UPstream::master())
        {
            for
            (
                int slave=UPstream::firstSlave();
                slave<=UPstream::lastSlave();
                slave++
            )
            {
                scalar value;

                if
                (
                    MPI_Recv
                    (
                        &value,
                        1,
                        MPI_SCALAR,
                        UPstream::procID(slave),
                        tag,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Recv failed"
                        << Foam::abort(FatalError);
                }

                Value = bop(Value, value);
            }
        }
        else
        {
            if
            (
                MPI_Send
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    UPstream::procID(UPstream::masterNo()),
                    tag,
                    MPI_COMM_WORLD
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Send failed"
                    << Foam::abort(FatalError);
            }
        }


        if (UPstream::master())
        {
            for
            (
                int slave=UPstream::firstSlave();
                slave<=UPstream::lastSlave();
                slave++
            )
            {
                if
                (
                    MPI_Send
                    (
                        &Value,
                        1,
                        MPI_SCALAR,
                        UPstream::procID(slave),
                        tag,
                        MPI_COMM_WORLD
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Send failed"
                        << Foam::abort(FatalError);
                }
            }
        }
        else
        {
            if
            (
                MPI_Recv
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    UPstream::procID(UPstream::masterNo()),
                    tag,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Recv failed"
                    << Foam::abort(FatalError);
            }
        }
    }
    else
    {
        scalar sum;
        MPI_Allreduce(&Value, &sum, 1, MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);
        Value = sum;

        /*
        int myProcNo = UPstream::myProcNo();
        int nProcs = UPstream::nProcs();

        //
        // receive from children
        //
        int level = 1;
        int thisLevelOffset = 2;
        int childLevelOffset = thisLevelOffset/2;
        int childProcId = 0;

        while
        (
            (childLevelOffset < nProcs)
         && (myProcNo % thisLevelOffset) == 0
        )
        {
            childProcId = myProcNo + childLevelOffset;

            scalar value;

            if (childProcId < nProcs)
            {
                if
                (
                    MPI_Recv
                    (
                        &value,
                        1,
                        MPI_SCALAR,
                        UPstream::procID(childProcId),
                        tag,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Recv failed"
                        << Foam::abort(FatalError);
                }

                Value = bop(Value, value);
            }

            level++;
            thisLevelOffset <<= 1;
            childLevelOffset = thisLevelOffset/2;
        }

        //
        // send and receive from parent
        //
        if (!UPstream::master())
        {
            int parentId = myProcNo - (myProcNo % thisLevelOffset);

            if
            (
                MPI_Send
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    UPstream::procID(parentId),
                    tag,
                    MPI_COMM_WORLD
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Send failed"
                    << Foam::abort(FatalError);
            }

            if
            (
                MPI_Recv
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    UPstream::procID(parentId),
                    tag,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Recv failed"
                    << Foam::abort(FatalError);
            }
        }


        //
        // distribute to my children
        //
        level--;
        thisLevelOffset >>= 1;
        childLevelOffset = thisLevelOffset/2;

        while (level > 0)
        {
            childProcId = myProcNo + childLevelOffset;

            if (childProcId < nProcs)
            {
                if
                (
                    MPI_Send
                    (
                        &Value,
                        1,
                        MPI_SCALAR,
                        UPstream::procID(childProcId),
                        tag,
                        MPI_COMM_WORLD
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Send failed"
                        << Foam::abort(FatalError);
                }
            }

            level--;
            thisLevelOffset >>= 1;
            childLevelOffset = thisLevelOffset/2;
        }
        */
    }

    if (Pstream::debug)
    {
        Pout<< "Foam::reduce : reduced value:" << Value << endl;
    }
}


Foam::label Foam::UPstream::nRequests()
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label i)
{
    if (i < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.setSize(i);
    }
}


void Foam::UPstream::waitRequests(const label start)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<MPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );

        if
        (
            MPI_Waitall
            (
                waitRequests.size(),
                waitRequests.begin(),
                MPI_STATUSES_IGNORE
            )
        )
        {
            FatalErrorIn
            (
                "UPstream::waitRequests()"
            )   << "MPI_Waitall returned with error" << Foam::endl;
        }

        resetRequests(start);
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorIn
        (
            "UPstream::finishedRequest(const label)"
        )   << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    int flag;
    MPI_Test
    (
       &PstreamGlobals::outstandingRequests_[i],
       &flag,
        MPI_STATUS_IGNORE
    );

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait for request:" << i
            << endl;
    }

    return flag != 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
