/*---------------------------------------------------------------------------*\
  L2-Roe low-dissipation Roe flux for low Mach numbers
  Implementation based on Oßwald et al. (2016) and Rieper (2011).

  Key corrections:
  - Local Mach number based on max(ML, MR), not Roe-averaged velocity
  - Scaling applied to velocity jumps in wave strength computation
  - Shock switch properly applied to both normal and tangential components
  - Physical flux evaluation unchanged
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  LMRoe - Low Mach number fix for Roe's approximate Riemann solver
  Implementation based on Rieper (2011) "A low-Mach number fix for Roe's 
  approximate Riemann solver", Journal of Computational Physics 230 (2011).

  Key differences from L2Roe:
  - Only the normal velocity jump is scaled by the local Mach number
  - Tangential velocity jumps remain UNSCALED
  - This is sufficient to fix the low Mach accuracy problem
  - L2Roe additionally scales tangential jumps to reduce dissipation
    at high wavenumbers (relevant for LES/DES)
\*---------------------------------------------------------------------------*/

#include "SGTNPRoeFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(SGTNPRoeFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, SGTNPRoeFlux, dictionary);
}

void Foam::SGTNPRoeFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf,
    const scalar& meshPhi,
    const tensor& R,
    const tensor& RTranspos,
    const scalar& fp1Left,
    const scalar& fp1Right
) const
{
    // cell face *normal* velocity w_n
    const scalar w_n = meshPhi / (magSf + VSMALL);


    // normal vector
    vector normalVector = Sf/magSf;


    //! Step 1: decode rho left and right:
    scalar rhoLeft  = (pLeft  +  pInf) / ((gamma - 1.0) * CvLeft * TLeft);
    scalar rhoRight = (pRight +  pInf) / ((gamma - 1.0) * CvRight * TRight);

    //! Thoenber modification of the reconstructed velocites
    // speeds of sound (left and right states)
    const scalar cLeft = Foam::sqrt(max((gamma*(pLeft + pInf))/rhoLeft,SMALL));
    const scalar cRight = Foam::sqrt(max((gamma*(pRight + pInf))/rhoRight,SMALL));
    
    // Local Mach numbers at left and right states (based on total velocity magnitude)
    const scalar ML = mag(ULeft  & normalVector) / (cLeft + VSMALL);
    const scalar MR = mag(URight & normalVector) / (cRight + VSMALL);
    
    // Scaling factor: min(1, max(ML, MR)) as per Equation (6) in the paper
    const scalar Ma_local = max(ML, MR); //!  increase

    // Info << "Ma_local= "<<Ma_local<<endl;
    const scalar zeta1 = min(1.0, Ma_local) ;


    // scalar fp = pow3(min( (pLeft/pRight) , (pRight/pLeft)));
    scalar fp = min(fp1Left,fp1Right);
    
    scalar zeta = 1- (1-zeta1) *fp ;

    const vector ULeft_star = 0.5*(ULeft + URight)   + zeta * 0.5*(ULeft - URight);
    const vector URight_star = 0.5*(ULeft + URight)  + zeta * 0.5*(URight - ULeft);


    const scalar rhoELeft   = rhoLeft*(((pLeft + gamma*pInf)/(pLeft + pInf) )*CvLeft*TLeft + q) + 0.5*rhoLeft*magSqr(ULeft_star);
    const scalar rhoERight  = rhoRight*(((pRight + gamma*pInf)/(pRight + pInf) )*CvRight*TRight + q) + 0.5*rhoRight*magSqr(URight_star);

        // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;


    // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
    // it must be computed for each cell and evaluated at each face
    // through reconstruction
    const scalar kappaLeft = gamma;
    const scalar kappaRight = gamma;


    // Compute left and right contravariant velocities:
    const scalar contrVLeft  = (ULeft_star & normalVector);
    const scalar contrVRight = (URight_star & normalVector);


    // // Step 2: compute Roe averged quantities for face:
    // const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1 - wLeft;

    const vector UTilde = ULeft_star*wLeft + URight_star*wRight;
    const scalar HTilde  = HLeft*wLeft + HRight*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

    // Roe Speed of sound
    const scalar cTilde =
        sqrt(max((gamma - 1)*(HTilde  - 0.5*qTildeSquare) -q, SMALL));

    // Roe averaged contravariant velocity
    const scalar contrVTilde = (UTilde & normalVector);

    // Step 3: compute primitive differences:
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = URight_star - ULeft_star;
    const scalar deltaContrV = (deltaU & normalVector);

    
    
    // // Apply scaling only away from shocks (ssw = 0)
    // const scalar zeta_eff = shockPresent ? 1.0 : zeta;
    // const scalar zeta_eff = zeta;
 
    // // Apply scaling only away from shocks (ssw = 0)
    // const scalar zeta_eff = shockPresent ? 1.0 : zeta;
    // const scalar zeta_eff = zeta;
 
        // -------------------------
    // Wave strengths using scaled jumps
    // -------------------------
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // // Acoustic wave strengths with SCALED normal velocity jump
    // const scalar r1 = (deltaP - rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    // const scalar r3 = (deltaP + rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    
    // // Entropy wave strength (unscaled)
    // const scalar r2 = deltaRho - deltaP/(cTilde*cTilde);

    const scalar r1 =
        (deltaP - rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));
    const scalar r2 = deltaRho - deltaP/sqr(cTilde);
    const scalar r3 =
        (deltaP + rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));

    // -------------------------
    // Eigenvectors (right eigenvectors of Roe matrix)
    // -------------------------
   // rho row:
    const scalar l1rho = 1;
    const scalar l2rho = 1;
    const scalar l3rho = 0;
    const scalar l4rho = 1;

    // first U column
    const vector l1U = UTilde - cTilde*normalVector;

    // second U column
    const vector l2U = UTilde;

    // third U column
    const vector l3U = deltaU - deltaContrV*normalVector;

    // fourth U column
    const vector l4U = UTilde + cTilde*normalVector;

    // E row
    const scalar l1e = HTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU) - contrVTilde*deltaContrV;
    const scalar l4e = HTilde + cTilde*contrVTilde;

    // Step 6: compute eigenvalues

    // -------------------------
    // Eigenvalues (wave speeds)
    // -------------------------
    
    scalar lambda1 = mag(contrVTilde - cTilde - w_n);
    scalar lambda2 = mag(contrVTilde - w_n);
    scalar lambda3 = mag(contrVTilde + cTilde - w_n);
    // Step 7: vanLeer entropy correction
    const scalar lambda1L = contrVTilde - cLeft  - w_n;
    const scalar lambda3L = contrVTilde + cLeft  - w_n;

    // Right state eigenvalues
    const scalar lambda1R = contrVTilde - cRight - w_n;
    const scalar lambda3R = contrVTilde + cRight - w_n;

    const scalar deltaLambda1 = max(lambda1R-lambda1L,0.0);
    const scalar deltaLabmda3 = max(lambda3R-lambda3L,0.0);
    
    if (lambda1 < 2*deltaLambda1)
    {
        lambda1 = (sqr(lambda1)/(4*deltaLambda1)) + deltaLambda1;

    }
    if (lambda3 < 2*deltaLabmda3)
    {
        lambda3 = (sqr(lambda3)/(4*deltaLabmda3)) + deltaLabmda3;

    }
    // Step 7: check for Harten entropy correction
    // -------------------------
    // Entropy correction (Harten-Hyman type)
    // -------------------------
 
    // const scalar UL = ULeft_star & normalVector;
    // const scalar UR = ULeft_star & normalVector;

    // // First eigenvalue: U - c
    // scalar eps = 2*max(0,(UR - cRight) - (UL - cLeft));
    // if (lambda1 < eps)
    // {
    //     lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
    // }

    // // Second eigenvalue: U
    // eps = 2*max(0, UR - UL);
    // if (lambda2 < eps)
    // {
    //     lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
    // }

    // // Third eigenvalue: U + c
    // eps = 2*max(0,(UR + cRight) - (UL + cLeft));
    // if (lambda3 < eps)
    // {
    //     lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
    // }

    // -------------------------
    // Flux difference components
    // -------------------------


    // Step 8: Compute flux differences

    // -------------------------
    // Flux difference components
    // -------------------------
    
    const scalar diffF11 = lambda1*r1*l1rho;
    const vector diffF124 = lambda1*r1*l1U;
    const scalar diffF15 = lambda1*r1*l1e;

    const scalar diffF21 = lambda2*(r2*l2rho + rhoTilde*l3rho);
    const vector diffF224 = lambda2*(r2*l2U + rhoTilde*l3U);
    const scalar diffF25 = lambda2*(r2*l2e + rhoTilde*l3e);

    const scalar diffF31 = lambda3*r3*l4rho;
    const vector diffF324 = lambda3*r3*l4U;
    const scalar diffF35 = lambda3*r3*l4e;
  
    // Step 9: compute left and right fluxes
    // -------------------------
    // Physical fluxes (unchanged from standard Roe)
    // -------------------------
    
    const scalar fluxLeft11 = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft_star*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = HLeft*fluxLeft11;
    // const scalar fluxLeft15 = (rhoELeft + pLeft)*contrVLeft;

    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight_star*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15 = HRight*fluxRight11;
    // const scalar fluxRight15 = (rhoERight + pRight)*contrVRight;

    // Step 10: compute face flux 5-vector
    // -------------------------
    // Face flux assembly (Roe flux with ALE correction)
    // -------------------------
    
    const scalar flux1 = 0.5*(fluxLeft11 + fluxRight11 
                              - (rhoLeft + rhoRight)*w_n 
                              - (diffF11 + diffF21 + diffF31));

    const vector flux24 = 0.5*(fluxLeft124 + fluxRight124 
                               - (rhoLeft*ULeft_star + rhoRight*URight_star)*w_n 
                               - (diffF124 + diffF224 + diffF324));

    const scalar flux5 = 0.5*(fluxLeft15 + fluxRight15 
                              - (rhoELeft + rhoERight )*w_n 
                              - (diffF15 + diffF25 + diffF35));

    // Scale by face area
    rhoFlux  = flux1*magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5*magSf;
}

// ************************************************************************* //
