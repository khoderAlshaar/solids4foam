    // Read field bounds
    dictionary fieldBounds = mesh().solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", p().dimensions(), SMALL);
    dimensionedScalar pMax("pMax", p().dimensions(), 1e10);
    fieldBounds.lookup(p().name()) >> pMin.value() >> pMax.value();

    // Temperature bounds
    dimensionedScalar TMin("TMin", T_.dimensions(), SMALL);
    dimensionedScalar TMax("TMax", T_.dimensions(), 1e10);
    fieldBounds.lookup(T_.name()) >> TMin.value() >> TMax.value();

    // Density bounds
    dimensionedScalar rhoMin("rhoMin", rho_.dimensions(), SMALL);
    dimensionedScalar rhoMax("rhoMax", rho_.dimensions(), 1e10);
    fieldBounds.lookup(rho_.name()) >> rhoMin.value() >> rhoMax.value();
