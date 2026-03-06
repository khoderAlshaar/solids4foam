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

#include "SGL2RoeFlux.H"
// #include "addToRunTimeSelectionTable.H"

// namespace Foam
// {
//     defineTypeNameAndDebug(SGL2RoeFlux, 0);
//     addToRunTimeSelectionTable(dbnsFlux, SGL2RoeFlux, dictionary);
// }

void Foam::SGL2RoeFlux::evaluateFlux
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
    const scalar& pInf,
    const scalar& q
) const
{

    const scalar gamma =  (RLeft / CvLeft) + 1 ;

    // cell face *normal* velocity w_n
    const scalar w_n = meshPhi / (magSf + VSMALL);

    //! Step 1: decode rho left and right:
    scalar rhoLeft  = (pLeft  +  pInf) / ((gamma - 1.0) * CvLeft * TLeft);
    scalar rhoRight = (pRight +  pInf) / ((gamma - 1.0) * CvRight * TRight);

    //! Decode left and right total energy:
    // total energy per unit mass
    // scalar eLeft  = (pLeft  + gamma*pInf)/((gamma-1)*rhoLeft)  + 0.5*magSqr(ULeft);
    // scalar eRight = (pRight + gamma*pInf)/((gamma-1)*rhoRight) + 0.5*magSqr(URight);

    // const scalar eLeft = (CvLeft*TLeft + pInf/rhoLeft +0.5*magSqr(ULeft));
    // const scalar eRight =(CvRight*TRight + pInf/rhoRight +0.5*magSqr(URight));

    
    // Compute left and right total enthalpies:
    // const scalar HLeft = eLeft + pLeft/rhoLeft;
    // const scalar HRight = eRight + pRight/rhoRight;

        // DensityTotalEnergy
    // const scalar rhoELeft = rhoLeft*(CvLeft*TLeft + pInf/rhoLeft +0.5*magSqr(ULeft));
    // const scalar rhoERight = rhoRight*(CvRight*TRight + pInf/rhoRight +0.5*magSqr(URight));  
    //       // DensityTotalEnergy
    // const scalar rhoELeft =  ((pLeft + gamma*pInf)/((gamma - 1)))+ rhoLeft*q  +0.5*rhoLeft*magSqr(ULeft);
    // const scalar rhoERight = ((pRight + gamma*pInf)/((gamma - 1)))+ rhoRight*q  +0.5*rhoRight*magSqr(URight);

    const scalar rhoELeft   = rhoLeft*(((pLeft + gamma*pInf)/(pLeft + pInf) )*CvLeft*TLeft + q) + 0.5*rhoLeft*magSqr(ULeft);
    const scalar rhoERight  = rhoRight*(((pRight + gamma*pInf)/(pRight + pInf) )*CvRight*TRight + q) + 0.5*rhoRight*magSqr(URight);

        // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;


    // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
    // it must be computed for each cell and evaluated at each face
    // through reconstruction
    const scalar kappaLeft = gamma;
    const scalar kappaRight = gamma;

    // normal vector
    vector normalVector = Sf/magSf;

    // Compute left and right contravariant velocities:
    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);


    // // Step 2: compute Roe averged quantities for face:
    // const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar HTilde  = HLeft*wLeft + HRight*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

    // Roe Speed of sound
    const scalar cTilde =
        sqrt(max((gamma - 1)*(HTilde  - 0.5*qTildeSquare), SMALL));

    // Roe averaged contravariant velocity
    const scalar contrVTilde = (UTilde & normalVector);

    // Step 3: compute primitive differences:
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = URight - ULeft;
    const scalar deltaContrV = (deltaU & normalVector);

    // Step 4: compute wave strengths:

    // -------------------------
    // L2-Roe scaling: Compute local Mach number based on LEFT and RIGHT states
    // -------------------------
    
    // speeds of sound (left and right states)
    // const scalar cLeft = Foam::sqrt(max((gamma*(pLeft + pInf))/rhoLeft,SMALL));
    // const scalar cRight = Foam::sqrt(max((gamma*(pRight + pInf))/rhoRight,SMALL));
    
    // // Local Mach numbers at left and right states (based on total velocity magnitude)
    // const scalar ML = mag(ULeft) / (cLeft + VSMALL);
    // const scalar MR = mag(URight) / (cRight + VSMALL);
    
    // // Scaling factor: min(1, max(ML, MR)) as per Equation (6) in the paper
    // const scalar Ma_local = max(ML, MR); //!  increase

    // // Info << "Ma_local= "<<Ma_local<<endl;
    // const scalar zeta = min(1.0, Ma_local) ;
    // -------------------------
    // Shock switch (Portela style with modifications)
    // -------------------------
    const vector UTilde_normal = contrVTilde * normalVector;
    const vector UTilde_tangent = UTilde - UTilde_normal;
    const scalar VTildeMag = mag(UTilde_tangent);
    
    // Local Mach number as per Rieper (2011), Eq. (3.16):
    // Ma_local = (|U_n| + |V_t|) / a
    const scalar Ma_local = (mag(UTilde_normal) + mag(UTilde_tangent)) / (cTilde + VSMALL);
    
    // Scaling factor: min(Ma_local, 1)
    const scalar zeta = min(1.0, Ma_local);

    // Simplified shock switch: active if significant pressure jump exists
    // Can be made more sophisticated based on pressure gradient
    // const scalar deltaP_threshold = 0.01*min(pLeft, pRight);
    // const bool shockPresent = (mag(deltaP) > deltaP_threshold);
    
    // // // Apply scaling only away from shocks (ssw = 0)
    // // const scalar zeta = shockPresent ? 1.0 : zeta;
    // const scalar zeta = zeta;

    // -------------------------
    // Scaled velocity jumps for L2-Roe
    // -------------------------
    
    // Scale normal velocity jump (LMRoe + L2Roe)
    const scalar deltaContrV_star = zeta * deltaContrV;
    
    // Tangential velocity components
    const vector deltaU_nVec = deltaContrV * normalVector;
    const vector deltaU_t = deltaU - deltaU_nVec;
    
    // Scale tangential velocity jump (L2Roe enhancement)
    const vector deltaU_t_star = zeta * deltaU_t;

        // -------------------------
    // Wave strengths using scaled jumps
    // -------------------------
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Acoustic wave strengths with SCALED normal velocity jump
    const scalar r1 = (deltaP - rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    const scalar r3 = (deltaP + rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    
    // Entropy wave strength (unscaled)
    const scalar r2 = deltaRho - deltaP/(cTilde*cTilde);


    // -------------------------
    // Eigenvectors (right eigenvectors of Roe matrix)
    // -------------------------
    
    const scalar l1rho = 1.0;
    const scalar l2rho = 1.0;
    const scalar l3rho = 0.0;
    const scalar l4rho = 1.0;

    const vector l1U = UTilde - cTilde*normalVector;
    const vector l2U = UTilde;
    const vector l3U = deltaU_t_star;  // Shear wave uses SCALED tangential jump
    const vector l4U = UTilde + cTilde*normalVector;

    const scalar l1e = HTilde  - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU_t_star);  // Consistent with scaled jump
    const scalar l4e = HTilde  + cTilde*contrVTilde;

    // Step 6: compute eigenvalues

    // -------------------------
    // Eigenvalues (wave speeds)
    // -------------------------
    
    scalar lambda1 = mag(contrVTilde - cTilde - w_n);
    scalar lambda2 = mag(contrVTilde - w_n);
    scalar lambda3 = mag(contrVTilde + cTilde - w_n);

    // Step 7: check for Harten entropy correction
    // -------------------------
    // Entropy correction (Harten-Hyman type)
    // -------------------------
    // Compute shock indicator
    // const scalar UL = contrVLeft;
    // const scalar UR = contrVRight;

    // scalar eps1 = 2.0*max(0.0, (UR - cRight) - (UL - cLeft));
    // scalar eps2 = 2.0*max(0.0, UR - UL);
    // scalar eps3 = 2.0*max(0.0, (UR + cRight) - (UL + cLeft));

    // if (lambda1 < eps1 && eps1 > VSMALL)
    // {
    //     lambda1 = (sqr(lambda1) + sqr(eps1))/(2.0*eps1);
    // }

    // if (lambda2 < eps2 && eps2 > VSMALL)
    // {
    //     lambda2 = (sqr(lambda2) + sqr(eps2))/(2.0*eps2);
    // }

    // if (lambda3 < eps3 && eps3 > VSMALL)
    // {
    //     lambda3 = (sqr(lambda3) + sqr(eps3))/(2.0*eps3);
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
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = HLeft*fluxLeft11;
    // const scalar fluxLeft15 = (rhoELeft + pLeft)*contrVLeft;

    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
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
                               - (rhoLeft*ULeft + rhoRight*URight)*w_n 
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
