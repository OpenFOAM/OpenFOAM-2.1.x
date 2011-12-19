dictionary additional = mesh.solutionDict().subDict("additional");

// pressure work term for enthalpy equation
bool pressureWorkTimeDerivative =
    additional.lookupOrDefault("pressureWorkTimeDerivative", false);

// flag to activate solve transport for each specie (Y vector)
bool solveSpecies = additional.lookupOrDefault("solveSpecies", true);
