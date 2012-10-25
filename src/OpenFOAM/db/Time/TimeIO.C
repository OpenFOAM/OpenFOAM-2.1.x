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

#include "Time.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Time::readDict()
{
    if (!deltaTchanged_)
    {
        deltaT_ = readScalar(controlDict_.lookup("deltaT"));
    }

    if (controlDict_.found("writeControl"))
    {
        writeControl_ = writeControlNames_.read
        (
            controlDict_.lookup("writeControl")
        );
    }

    scalar oldWriteInterval = writeInterval_;
    scalar oldSecondaryWriteInterval = secondaryWriteInterval_;

    if (controlDict_.readIfPresent("writeInterval", writeInterval_))
    {
        if (writeControl_ == wcTimeStep && label(writeInterval_) < 1)
        {
            FatalIOErrorIn("Time::readDict()", controlDict_)
                << "writeInterval < 1 for writeControl timeStep"
                << exit(FatalIOError);
        }
    }
    else
    {
        controlDict_.lookup("writeFrequency") >> writeInterval_;
    }


    // Additional writing
    if (controlDict_.found("secondaryWriteControl"))
    {
        secondaryWriteControl_ = writeControlNames_.read
        (
            controlDict_.lookup("secondaryWriteControl")
        );

        if
        (
            controlDict_.readIfPresent
            (
                "secondaryWriteInterval",
                secondaryWriteInterval_
            )
        )
        {
            if
            (
                secondaryWriteControl_ == wcTimeStep
             && label(secondaryWriteInterval_) < 1
            )
            {
                FatalIOErrorIn("Time::readDict()", controlDict_)
                    << "secondaryWriteInterval < 1"
                    << " for secondaryWriteControl timeStep"
                    << exit(FatalIOError);
            }
        }
        else
        {
            controlDict_.lookup("secondaryWriteFrequency")
                >> secondaryWriteInterval_;
        }
    }



    if (oldWriteInterval != writeInterval_)
    {
        switch (writeControl_)
        {
            case wcRunTime:
            case wcAdjustableRunTime:
                // Recalculate outputTimeIndex_ to be in units of current
                // writeInterval.
                outputTimeIndex_ = label
                (
                    outputTimeIndex_
                  * oldWriteInterval
                  / writeInterval_
                );
            break;

            default:
            break;
        }
    }
    if (oldSecondaryWriteInterval != secondaryWriteInterval_)
    {
        switch (secondaryWriteControl_)
        {
            case wcRunTime:
            case wcAdjustableRunTime:
                // Recalculate secondaryOutputTimeIndex_ to be in units of
                // current writeInterval.
                secondaryOutputTimeIndex_ = label
                (
                    secondaryOutputTimeIndex_
                  * oldSecondaryWriteInterval
                  / secondaryWriteInterval_
                );
            break;

            default:
            break;
        }
    }

    if (controlDict_.readIfPresent("purgeWrite", purgeWrite_))
    {
        if (purgeWrite_ < 0)
        {
            WarningIn("Time::readDict()")
                << "invalid value for purgeWrite " << purgeWrite_
                << ", should be >= 0, setting to 0"
                << endl;

            purgeWrite_ = 0;
        }
    }

    if (controlDict_.found("timeFormat"))
    {
        const word formatName(controlDict_.lookup("timeFormat"));

        if (formatName == "general")
        {
            format_ = general;
        }
        else if (formatName == "fixed")
        {
            format_ = fixed;
        }
        else if (formatName == "scientific")
        {
            format_ = scientific;
        }
        else
        {
            WarningIn("Time::readDict()")
                << "unsupported time format " << formatName
                << endl;
        }
    }

    controlDict_.readIfPresent("timePrecision", precision_);

    // stopAt at 'endTime' or a specified value
    // if nothing is specified, the endTime is zero
    if (controlDict_.found("stopAt"))
    {
        stopAt_ = stopAtControlNames_.read(controlDict_.lookup("stopAt"));

        if (stopAt_ == saEndTime)
        {
            controlDict_.lookup("endTime") >> endTime_;
        }
        else
        {
            endTime_ = GREAT;
        }
    }
    else if (!controlDict_.readIfPresent("endTime", endTime_))
    {
        endTime_ = 0;
    }

    dimensionedScalar::name() = timeName(value());

    if (controlDict_.found("writeVersion"))
    {
        writeVersion_ = IOstream::versionNumber
        (
            controlDict_.lookup("writeVersion")
        );
    }

    if (controlDict_.found("writeFormat"))
    {
        writeFormat_ = IOstream::formatEnum
        (
            controlDict_.lookup("writeFormat")
        );
    }

    if (controlDict_.found("writePrecision"))
    {
        IOstream::defaultPrecision
        (
            readUint(controlDict_.lookup("writePrecision"))
        );

        Sout.precision(IOstream::defaultPrecision());
        Serr.precision(IOstream::defaultPrecision());

        Pout.precision(IOstream::defaultPrecision());
        Perr.precision(IOstream::defaultPrecision());

        FatalError().precision(IOstream::defaultPrecision());
        FatalIOError.error::operator()().precision
        (
            IOstream::defaultPrecision()
        );
    }

    if (controlDict_.found("writeCompression"))
    {
        writeCompression_ = IOstream::compressionEnum
        (
            controlDict_.lookup("writeCompression")
        );
    }

    controlDict_.readIfPresent("graphFormat", graphFormat_);
    controlDict_.readIfPresent("runTimeModifiable", runTimeModifiable_);

    if (!runTimeModifiable_ && controlDict_.watchIndex() != -1)
    {
        removeWatch(controlDict_.watchIndex());
        controlDict_.watchIndex() = -1;
    }
}


bool Foam::Time::read()
{
    if (controlDict_.regIOobject::read())
    {
        readDict();
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::Time::readModifiedObjects()
{
    if (runTimeModifiable_)
    {
        // Get state of all monitored objects (=registered objects with a
        // valid filePath).
        // Note: requires same ordering in objectRegistries on different
        // processors!
        monitorPtr_().updateStates
        (
            (
                regIOobject::fileModificationChecking == inotifyMaster
             || regIOobject::fileModificationChecking == timeStampMaster
            ),
            Pstream::parRun()
        );

        // Time handling is special since controlDict_ is the one dictionary
        // that is not registered to any database.

        if (controlDict_.readIfModified())
        {
           readDict();
           functionObjects_.read();
        }

        bool registryModified = objectRegistry::modified();

        if (registryModified)
        {
            objectRegistry::readModifiedObjects();
        }
    }
}


bool Foam::Time::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    if (outputTime())
    {
        const word tmName(timeName());

        IOdictionary timeDict
        (
            IOobject
            (
                "time",
                tmName,
                "uniform",
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        timeDict.add("value", value());
        timeDict.add("name", string(tmName));
        timeDict.add("index", timeIndex_);
        timeDict.add("deltaT", deltaT_);
        timeDict.add("deltaT0", deltaT0_);

        timeDict.regIOobject::writeObject(fmt, ver, cmp);
        bool writeOK = objectRegistry::writeObject(fmt, ver, cmp);

        if (writeOK && purgeWrite_)
        {
            previousOutputTimes_.push(tmName);

            while (previousOutputTimes_.size() > purgeWrite_)
            {
                rmDir(objectRegistry::path(previousOutputTimes_.pop()));
            }
        }

        return writeOK;
    }
    else
    {
        return false;
    }
}


bool Foam::Time::writeNow()
{
    outputTime_ = true;
    return write();
}


bool Foam::Time::writeAndEnd()
{
    stopAt_  = saWriteNow;
    endTime_ = value();

    return writeNow();
}


void Foam::Time::writeOnce()
{
    writeOnce_ = true;
}


// ************************************************************************* //
