The 0/ field files contain nonsense patchFields. All interesting
work is done using the changeDictionaryDicts.



Note:To run with directional thermo:

- compile chtMultiRegionFoam/solid/setRegionSolidFields.H with

    tmp<volSymmTensorField> tkappa = thermo.directionalK();
    const volSymmTensorField& kappa = tkappa();

- change in e.g. heater:

    - in constant/heater/solidThermophysicalProperties:

        thermoType directionalKSolidThermo;

    - in 0/heater/T:

        K               directionalSolidThermo;
