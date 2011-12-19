dictionary additional = mesh.solutionDict().subDict("additional");

// pressure work term for enthalpy equation
bool pressureWork = additional.lookupOrDefault("pressureWork", true);
bool pressureWorkTimeDerivative =
    additional.lookupOrDefault("pressureWorkTimeDerivative", true);

// flag to activate solve transport for each specie (Y vector)
bool solveSpecies = additional.lookupOrDefault("solveSpecies", true);
