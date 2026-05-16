#include "ConvergenceControl.H"
#include "utils.H"
#include "volFields.H"
#include "OFstream.H"
#include "Pstream.H"
#include "OSspecific.H"

// =======================
// Constructor
// =======================

ConvergenceControl::ConvergenceControl
(
    const Time& runTime,
    const fvMesh& mesh
)
:
    runTime_(runTime),
    mesh_(mesh),

    pseudoMaxIters_(50),
    pseudoAbsTol_(1e-6),
    pseudoRelTol_(0.1),
    pseudoStagnationTol_(1e-8),
    pseudoStagnationLimit_(5),
    infoFrequency_(10),

    initialPseudoResidual_(GREAT),
    previousPseudoResidual_(GREAT),
    pseudoStagnationCount_(0),

    physicalAbsTol_(1e-6),
    physicalConsecutiveTol_(5),
    physicalTolCounter_(0),

    solverType_("default"),
    pseudoFilePtr_(),
    physicalFilePtr_()
{
    if (Pstream::master())
    {
        Info<< "Creating residual files" << endl;

        fileName outDir;

        if (Pstream::parRun())
        {
            outDir = runTime_.path()/"..";
        }
        else
        {
            outDir = runTime_.path();
        }

        outDir = outDir/"postProcessing"/runTime_.timeName();
        mkDir(outDir);

        // --- Pseudo residuals file
        pseudoFilePtr_.set
        (
            new OFstream(outDir/"pseudoResiduals.dat")
        );

        pseudoFilePtr_()
            << "# timeIndex time pseudoIter residual rhoRes rhoURes rhoERes\n";

        // --- Physical residuals file
        physicalFilePtr_.set
        (
            new OFstream(outDir/"physicalResiduals.dat")
        );

        physicalFilePtr_()
            << "# timeIndex time rhoRes rhoURes rhoERes\n";
    }
}

// =======================
// Dictionary read
// =======================

void ConvergenceControl::read(const dictionary& dict)
{
    if (dict.found("multiStage"))
    {
        const dictionary& ms = dict.subDict("multiStage");

        ms.lookup("pseudoMaxIters") >> pseudoMaxIters_;
        ms.lookup("pseudoAbsTol")   >> pseudoAbsTol_;

        infoFrequency_ = lookupOrDefault<int>(ms, "infoFrequency", infoFrequency_);

        pseudoStagnationTol_ =
            lookupOrDefault<scalar>(ms, "pseudoStagnationTol", pseudoStagnationTol_);

        pseudoStagnationLimit_ =
            lookupOrDefault<label>(ms, "pseudoStagnationLimit", pseudoStagnationLimit_);

        ms.lookup("solverType") >> solverType_;
    }

    if (dict.found("globalConvergence"))
    {
        const dictionary& gc = dict.subDict("globalConvergence");

        physicalAbsTol_ =
            lookupOrDefault<scalar>(gc, "physicalAbsTol", physicalAbsTol_);

        physicalConsecutiveTol_ =
            lookupOrDefault<label>(gc, "physicalConsecutiveTol", physicalConsecutiveTol_);
    }
}

// =======================
// Append functions
// =======================

void ConvergenceControl::appendPseudoResidual
(
    int iter,
    scalar res,
    scalar rhoRes,
    scalar rhoURes,
    scalar rhoERes
)
{
    if (pseudoFilePtr_.valid())
    {
        pseudoFilePtr_()
            << runTime_.timeIndex() << " "
            << runTime_.value()     << " "
            << iter                 << " "
            << res                  << " "
            << rhoRes               << " "
            << rhoURes              << " "
            << rhoERes              << "\n";

        pseudoFilePtr_().flush();
    }
}


void ConvergenceControl::appendPhysicalResidual
(
    scalar dR,
    scalar dRU,
    scalar dRE
)
{
    if (physicalFilePtr_.valid())
    {
        physicalFilePtr_()
            << runTime_.timeIndex() << " "
            << runTime_.value()     << " "
            << dR                   << " "
            << dRU                  << " "
            << dRE                  << "\n";

        physicalFilePtr_().flush();
    }
}


// =======================
// Residual definitions
// =======================
template<class FieldType>
static inline scalar computePseudoResidualT(const FieldType& f)
{
    // denom based on max change from previous physical time (oldTime)
    scalar denom = gMax(mag(f.internalField() - f.oldTime().internalField()));
    if (denom < SMALL)
    {
        denom = max(gMax(mag(f.internalField())), SMALL);
    }

    return gMax(mag(f.internalField() - f.prevIter().internalField())) / denom;
}



// template<class FieldType>
// static inline scalar computePseudoResidualT(const FieldType& f)
// {
//     const FieldType& fPrev = f.prevIter();

//     return Foam::sqrt
//     (
//         gAverage(magSqr(f.internalField() - fPrev.internalField()))
//     );
// }


template<class FieldType>
static inline scalar computePhysicalResidualT(const FieldType& f)
{
    const FieldType& f0 = f.oldTime();

    return Foam::sqrt
    (
        gAverage(magSqr(f.internalField() - f0.internalField()))
    );
}


// =======================
// PSEUDO convergence
// =======================

bool ConvergenceControl::pseudoConverged
(
    int iter,
    const dimensionedScalar& dtPseudo,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE
)
{
    scalar rhoRes  = computePseudoResidualT(rho);
    scalar rhoURes = computePseudoResidualT(rhoU);
    scalar rhoERes = computePseudoResidualT(rhoE);

    // static scalar rho0 = -1, rhoU0 = -1, rhoE0 = -1;

    // if (iter == 1)
    // {
    //     rho0  = max(rhoRes, SMALL);
    //     rhoU0 = max(rhoURes, SMALL);
    //     rhoE0 = max(rhoERes, SMALL);
    // }

    // rhoRes  /= rho0;
    // rhoURes /= rhoU0;
    // rhoERes /= rhoE0;

    scalar residual = max(max(rhoRes, rhoURes), rhoERes);

    if (iter == 1)
    {
        initialPseudoResidual_ = residual;

        Info<< "    Iter  Res(max)   rhoRes   rhoURes   rhoERes   dtPseudo\n";
    }

    previousPseudoResidual_ = residual;

    // --- Write every iteration (important!)
    // appendPseudoResidual(iter, residual, rhoRes, rhoURes, rhoERes);

    bool conv = false;

    if (residual < pseudoAbsTol_)
    {
        Info<< "Pseudo converged (abs tol)\n";
        conv = true;
    }

    if (!conv && iter >= pseudoMaxIters_)
    {
        Warning<< "Pseudo max iterations reached\n";
        conv = true;
    }

    if (iter % infoFrequency_ == 0 || conv)
    {
        Info<< "    " << iter
            << "  " << residual
            << "  " << rhoRes
            << "  " << rhoURes
            << "  " << rhoERes
            << "  " << dtPseudo.value()
            << endl;
    }

    return conv;
}


// =======================
// PHYSICAL convergence
// =======================

bool ConvergenceControl::physicalConverged
(
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE
)
{
    scalar L2rho  = computePhysicalResidualT(rho);
    scalar L2rhoU = computePhysicalResidualT(rhoU);
    scalar L2rhoE = computePhysicalResidualT(rhoE);

    // static scalar rho0 = -1, rhoU0 = -1, rhoE0 = -1;

    // if (runTime_.timeIndex() == 1) // first physical step
    // {
    //     rho0  = max(L2rho, SMALL);
    //     rhoU0 = max(L2rhoU, SMALL);
    //     rhoE0 = max(L2rhoE, SMALL);
    // }

    //  L2rho  = L2rho  / rho0;
    //  L2rhoU = L2rhoU / rhoU0;
    //  L2rhoE = L2rhoE / rhoE0;


    scalar res = max(max(L2rho, L2rhoU), L2rhoE);

    appendPhysicalResidual(L2rho, L2rhoU, L2rhoE);

    if (res < physicalAbsTol_)
    {
        physicalTolCounter_++;
    }
    else
    {
        physicalTolCounter_ = 0;
    }

    if (physicalTolCounter_ >= physicalConsecutiveTol_)
    {
        Info<< "Global steady-state reached.\n";
        return true;
    }

    return false;
}