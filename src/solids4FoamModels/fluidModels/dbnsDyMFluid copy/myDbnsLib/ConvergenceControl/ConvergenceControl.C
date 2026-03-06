#include "ConvergenceControl.H"
#include "utils.H"     // your helper utilities, kept as before
#include "volFields.H" // ensure vol fields ops are available

// Constructor: open files and write header
ConvergenceControl::ConvergenceControl(const Time& runTime, const fvMesh& mesh )
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

    pseudoFile_(runTime_.path()/"pseudoResiduals.csv"),
    physicalFile_(runTime_.path()/"physicalResiduals.csv") // cleaned filename
{
    // Headers
    pseudoFile_ << "timeIndex,pseudoIter,residual\n";
    pseudoFile_.flush();

    physicalFile_ << "timeIndex,rho,rhoU,rhoE\n";
    physicalFile_.flush();
}


void ConvergenceControl::read(const dictionary& dict)
{
    if (dict.found("multiStage"))
    {
        const dictionary& ms = dict.subDict("multiStage");

        pseudoMaxIters_ = lookupOrDefault<label>(ms, "pseudoMaxIters", pseudoMaxIters_);
        pseudoAbsTol_   = lookupOrDefault<scalar>(ms, "pseudoAbsTol",   pseudoAbsTol_);
        pseudoRelTol_   = lookupOrDefault<scalar>(ms, "pseudoRelTol",   pseudoRelTol_);

        // alternative names / backwards compat
        pseudoMaxIters_ = lookupOrDefault<label>(ms, "numberSubCycles", pseudoMaxIters_);
        pseudoAbsTol_   = lookupOrDefault<scalar>(ms, "tolerance",      pseudoAbsTol_);
        pseudoRelTol_   = lookupOrDefault<scalar>(ms, "relTol",         pseudoRelTol_);
        infoFrequency_  = lookupOrDefault<int>(ms, "infoFrequency", infoFrequency_);

        pseudoStagnationTol_   = lookupOrDefault<scalar>(ms, "pseudoStagnationTol",   pseudoStagnationTol_);
        pseudoStagnationLimit_ = lookupOrDefault<label>(ms,  "pseudoStagnationLimit", pseudoStagnationLimit_);
    }

    if (dict.found("globalConvergence"))
    {
        const dictionary& gc = dict.subDict("globalConvergence");

        physicalAbsTol_         = lookupOrDefault<scalar>(gc, "physicalAbsTol", physicalAbsTol_);
        physicalConsecutiveTol_ = lookupOrDefault<label>(gc, "physicalConsecutiveTol", physicalConsecutiveTol_);
    }
}


// =======================
// On-the-fly appenders
// =======================

void ConvergenceControl::appendPseudoResidual(int iter, scalar res)
{
    pseudoFile_ << runTime_.timeIndex() << "," << iter << "," << res << "\n";
    pseudoFile_.flush();
}

void ConvergenceControl::appendPhysicalResidual(scalar dR, scalar dRU, scalar dRE)
{
    // Use the actual physical time index, not the counter of consecutive converged steps.
    physicalFile_ << runTime_.timeIndex() << "," << runTime_.value() << "," << dR << "," << dRU << "," << dRE << "\n";
    physicalFile_.flush();
}  


// -----------------------
// Templated residual helper (internal to this translation unit)
// Works for volScalarField and volVectorField (and other field types
// for which mag(), internalField(), oldTime(), prevIter() are defined).
// -----------------------
template<class FieldType>
static inline scalar computeResidualT(const FieldType& f)
{
    // denom based on max change from previous physical time (oldTime)
    scalar denom = gMax(mag(f.internalField() - f.oldTime().internalField()));
    if (denom < SMALL)
    {
        denom = max(gMax(mag(f.internalField())), SMALL);
    }

    return gMax(mag(f.internalField() - f.prevIter().internalField())) / denom;
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
    bool conv = false;

    // --- compute residuals using the templated helper
    scalar rhoRes  = computeResidualT(rho);
    scalar rhoURes = computeResidualT(rhoU);
    scalar rhoERes = computeResidualT(rhoE);

    // combined residual (use max to ensure all are converged)
    scalar residual = max(max(rhoRes, rhoURes), rhoERes);

    // initial residual
    if (iter == 0)
        initialPseudoResidual_ = residual;

    // stagnation
    if (previousPseudoResidual_ < GREAT)
    {
        scalar delta = mag(residual - previousPseudoResidual_);
        if (delta < pseudoStagnationTol_)
        {
            pseudoStagnationCount_++;
            if (pseudoStagnationCount_ >= pseudoStagnationLimit_)
            {
                Info<< "Pseudo converged by stagnation\n";
                conv = true;
                appendPseudoResidual(iter, residual);
            }
        }
        else
        {
            pseudoStagnationCount_ = 0;
        }
    }

    previousPseudoResidual_ = residual;

    // absolute tolerance
    if (!conv && residual < pseudoAbsTol_)
    {
        Info<< "Pseudo converged by absolute tolerance\n";
        conv = true;
        appendPseudoResidual(iter, residual);
    }

    // max iterations
    if (!conv && iter >= pseudoMaxIters_)
    {
        Warning<< "Pseudo max iterations reached\n";
        conv = true;
        appendPseudoResidual(iter, residual);
    }

    // logging
    if (iter == 0)
    {
        Info<< "    Corr, res(max), rhoRes, rhoURes, rhoERes, dtPseudo\n";
    }

    if (iter % infoFrequency_ == 0 || conv)
    {
        Info<< "    " << iter
            << ", " << residual
            << ", " << rhoRes
            << ", " << rhoURes
            << ", " << rhoERes
            << ", " << dtPseudo.value() << endl;
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
    const volScalarField& rho0 = rho.oldTime();

    scalar L2rho =
        Foam::sqrt
        (
            gSum
            (
                magSqr(rho.internalField() - rho0.internalField())
            * mesh_.V()  // multiply by cell volumes
            )
        );

    const volVectorField& rhoU0 = rhoU.oldTime();

    scalar L2rhoU =
        Foam::sqrt
        (
            gSum
            (
                magSqr(rhoU.internalField() - rhoU0.internalField())
            * mesh_.V()
            )
        );

    const volScalarField& rhoE0 = rhoE.oldTime();

    scalar L2rhoE =
        Foam::sqrt
        (
            gSum
            (
                magSqr(rhoE.internalField() - rhoE0.internalField())
            * mesh_.V()
            )
        );

const scalar eps = SMALL; // or 1e-16

scalar denomRho  = Foam::sqrt(gSum(magSqr(rho0.internalField())  * mesh_.V()) + eps);
scalar denomRhoU = Foam::sqrt(gSum(magSqr(rhoU0.internalField()) * mesh_.V()) + eps);
scalar denomRhoE = Foam::sqrt(gSum(magSqr(rhoE0.internalField()) * mesh_.V()) + eps);

scalar relL2rho  = L2rho  / denomRho;
scalar relL2rhoU = L2rhoU / denomRhoU;
scalar relL2rhoE = L2rhoE / denomRhoE;

    // choose combined metric (here L2 max among variables)
    scalar res = max(max(L2rho, L2rhoU), L2rhoE);

    // append relative L2 values to file for monitoring
    appendPhysicalResidual(relL2rho, relL2rhoU, relL2rhoE);

    // Info<< "Physical residual -> density = " << L2rho << endl;
    // Info<< "Physical residual -> linear momentum = " << L2rhoU << endl;
    // Info<< "Physical residual -> total energy = " << L2rhoE << endl;

    if (res < physicalAbsTol_)
        physicalTolCounter_++;
    else
        physicalTolCounter_ = 0;

    if (physicalTolCounter_ >= physicalConsecutiveTol_)
    {
        Info<< "Global steady-state reached.\n";
        return true;
    }

    return false;
}
