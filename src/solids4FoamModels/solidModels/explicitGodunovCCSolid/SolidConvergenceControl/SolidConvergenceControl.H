#ifndef CONVERGENCECONTROL_H
#define CONVERGENCECONTROL_H

#include "fvCFD.H"
#include "pointFields.H"
#include "pointMesh.H"

class SolidConvergenceControl
{
private:

    const Time& runTime_;
    const fvMesh& mesh_;

    // -------- pseudo-time parameters --------
    label  pseudoMaxIters_;
    scalar pseudoAbsTol_;
    scalar pseudoRelTol_;
    scalar pseudoStagnationTol_;
    label  pseudoStagnationLimit_;
    int    infoFrequency_;

    // internal tracking
    scalar initialPseudoResidual_;
    scalar previousPseudoResidual_;
    label  pseudoStagnationCount_;

    // -------- physical/global parameters --------
    scalar physicalAbsTol_;
    label  physicalConsecutiveTol_;

    // counter
    label physicalTolCounter_;

    // -------- Residual output files --------
    autoPtr<OFstream> pseudoFilePtr_;
    autoPtr<OFstream> physicalFilePtr_;

    // Solver type (MSMS, MSSS, SSSS)
    word solverType_;

public:

    // Constructor
    SolidConvergenceControl(const Time& runTime, const fvMesh& mesh);

    // Read from dictionary
    void read(const dictionary& dict);

    // =======================
    // Convergence checks
    // =======================
    bool pseudoConverged
    (
        int iter,
        const dimensionedScalar& dtPseudo,
        const volVectorField& lm,
        const volTensorField& P,
        const pointVectorField& xN
    );

    bool physicalConverged
    (
        const volVectorField& lm,
        const volTensorField& P,
        const pointVectorField& xN
    );

    // =======================
    // Residual writers
    // =======================

    // Pseudo residual (inner iterations)
    void appendPseudoResidual
    (
        int iter,
        scalar res,
        scalar rhoRes,
        scalar rhoURes,
        scalar rhoERes
    );

    // Physical residual (outer iterations / time)
    void appendPhysicalResidual
    (
        scalar dR,
        scalar dRU,
        scalar dRE
    );

    // =======================
    // Getters
    // =======================

    inline label& pseudoMaxIters()
    {
        return pseudoMaxIters_;
    }

    inline word& solverType()
    {
        return solverType_;
    }
};

#endif