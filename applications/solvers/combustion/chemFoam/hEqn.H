{
    if (constProp == "volume")
    {
        hs[0] = u0 + p[0]/rho[0] + integratedHeat;
    }
    else
    {
        hs[0] = hs0 + integratedHeat;
    }
}