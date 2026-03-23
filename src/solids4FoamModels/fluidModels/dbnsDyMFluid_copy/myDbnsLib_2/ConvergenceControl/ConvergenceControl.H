#ifndef CONVERGENCECONTROL_H
#define CONVERGENCECONTROL_H

#include "fvCFD.H"

class ConvergenceControl
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
    int infoFrequency_;

    // internal tracking
    scalar initialPseudoResidual_;
    scalar previousPseudoResidual_;
    label  pseudoStagnationCount_;

    // -------- physical/global parameters --------
    scalar physicalAbsTol_;
    label  physicalConsecutiveTol_;

    // counter
    label physicalTolCounter_;

    // -------- On-the-fly residual log files --------
    OFstream pseudoFile_;
    OFstream physicalFile_;

        word solverType_;   

public:

    // Constructor
    ConvergenceControl(const Time& runTime, const fvMesh& mesh );

    // Read from dictionary
    void read(const dictionary& dict);

    // Convergence checks
    bool pseudoConverged(int iter,
                         const dimensionedScalar& dtPseudo,
                         const volScalarField& rho,
                         const volVectorField& rhoU,
                         const volScalarField& rhoE);

    bool physicalConverged(const volScalarField& rho,
                           const volVectorField& rhoU,
                           const volScalarField& rhoE);

    // appenders
    void appendPseudoResidual(int iter, scalar res);
    void appendPhysicalResidual(scalar dR, scalar dRU, scalar dRE);

    // getter
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
