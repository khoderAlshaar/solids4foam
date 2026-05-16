/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hlleTNPFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(hlleTNPFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, hlleTNPFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hlleTNPFlux::evaluateFlux 
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
    // --------------------------------------------------
    // 1. Geometry
    // --------------------------------------------------
    const vector n = Sf / magSf;

    // --------------------------------------------------
    // 2. Thermodynamics
    // --------------------------------------------------
    const scalar kappaLeft  = (RLeft  + CvLeft)  / CvLeft;
    const scalar kappaRight = (RRight + CvRight) / CvRight;

    const scalar rhoLeft  = pLeft  / (RLeft  * TLeft);
    const scalar rhoRight = pRight / (RRight * TRight);

    const scalar aLeft  = Foam::sqrt(max(0.0, kappaLeft  * pLeft  / rhoLeft));
    const scalar aRight = Foam::sqrt(max(0.0, kappaRight * pRight / rhoRight));

    // --------------------------------------------------
    // 3. TNP velocity reconstruction (from paper)
    // --------------------------------------------------
    const scalar qL = (ULeft  & n);
    const scalar qR = (URight & n);

    const scalar ML = mag(qL) / (aLeft  + VSMALL);
    const scalar MR = mag(qR) / (aRight + VSMALL);

    const scalar zeta1 = min(1.0, max(ML, MR));
    const scalar fp = min(fp1Left, fp1Right);
    const scalar zeta = 1.0 - (1.0 - zeta1) * fp;

    const vector ULeft_star  = 0.5*(ULeft + URight) + 0.5*zeta*(ULeft - URight);
    const vector URight_star = 0.5*(ULeft + URight) + 0.5*zeta*(URight - ULeft);

    // --------------------------------------------------
    // 4. Conservative variables (reconstructed)
    // --------------------------------------------------
    const vector rhoULeft  = rhoLeft  * ULeft_star;
    const vector rhoURight = rhoRight * URight_star;

    const scalar rhoELeft =
        rhoLeft * (CvLeft * TLeft + 0.5 * magSqr(ULeft_star));

    const scalar rhoERight =
        rhoRight * (CvRight * TRight + 0.5 * magSqr(URight_star));

    const scalar qLeft  = ULeft_star  & n;
    const scalar qRight = URight_star & n;

    // --------------------------------------------------
    // 5. Roe averages (for wave speeds)
    // --------------------------------------------------
    const scalar rL = Foam::sqrt(max(0.0, rhoLeft));
    const scalar rR = Foam::sqrt(max(0.0, rhoRight));

    const scalar wL = rL / stabilise(rL + rR, VSMALL);
    const scalar wR = 1.0 - wL;

    const vector UTilde = wL * ULeft_star + wR * URight_star;
    const scalar unTilde = UTilde & n;

    const scalar HLeft  = (rhoELeft  + pLeft)  / rhoLeft;
    const scalar HRight = (rhoERight + pRight) / rhoRight;

    const scalar HTilde = wL * HLeft + wR * HRight;
    const scalar kappaTilde = wL * kappaLeft + wR * kappaRight;

    const scalar aTilde = Foam::sqrt
    (
        max(0.0, (kappaTilde - 1.0) * (HTilde - 0.5 * magSqr(UTilde)))
    );

    // --------------------------------------------------
    // 6. HLLE wave speeds (Einfeldt)
    // --------------------------------------------------
    const scalar SLeft  = min(qLeft - aLeft,  unTilde - aTilde);
    const scalar SRight = max(qRight + aRight, unTilde + aTilde);

    // --------------------------------------------------
    // 7. HLLE flux (NO star region!)
    // --------------------------------------------------
    const scalar invDen = 1.0 / stabilise(SRight - SLeft, VSMALL);

    // Left fluxes
    const scalar FL_rho = rhoLeft * qLeft;
    const vector FL_rhoU = rhoULeft * qLeft + pLeft * n;
    const scalar FL_rhoE = (rhoELeft + pLeft) * qLeft;

    // Right fluxes
    const scalar FR_rho = rhoRight * qRight;
    const vector FR_rhoU = rhoURight * qRight + pRight * n;
    const scalar FR_rhoE = (rhoERight + pRight) * qRight;

    // HLLE formula
    rhoFlux =
    (
        SRight * FL_rho
      - SLeft  * FR_rho
      + SRight * SLeft * (rhoRight - rhoLeft)
    ) * invDen * magSf;

    rhoUFlux =
    (
        SRight * FL_rhoU
      - SLeft  * FR_rhoU
      + SRight * SLeft * (rhoURight - rhoULeft)
    ) * invDen * magSf;

    rhoEFlux =
    (
        SRight * FL_rhoE
      - SLeft  * FR_rhoE
      + SRight * SLeft * (rhoERight - rhoELeft)
    ) * invDen * magSf;
}

// ************************************************************************* //

void Foam::hlleTNPFlux::evaluateFlux_Wall
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
    // --------------------------------------------------
    // 1. Geometry
    // --------------------------------------------------
    const vector n = Sf / magSf;


    
    const vector URight_wall = ULeft - 2.0*(ULeft & n)*n;


    // --------------------------------------------------
    // 2. Thermodynamics
    // --------------------------------------------------
    const scalar kappaLeft  = (RLeft  + CvLeft)  / CvLeft;
    const scalar kappaRight = (RRight + CvRight) / CvRight;

    const scalar rhoLeft  = pLeft  / (RLeft  * TLeft);
    const scalar rhoRight = pRight / (RRight * TRight);

    const scalar aLeft  = Foam::sqrt(max(0.0, kappaLeft  * pLeft  / rhoLeft));
    const scalar aRight = Foam::sqrt(max(0.0, kappaRight * pRight / rhoRight));

    // --------------------------------------------------
    // 3. TNP velocity reconstruction (from paper)
    // --------------------------------------------------
    const scalar qL = (ULeft  & n);
    const scalar qR = (URight_wall & n);

    const scalar ML = mag(qL) / (aLeft  + VSMALL);
    const scalar MR = mag(qR) / (aRight + VSMALL);

    const scalar zeta1 = min(1.0, max(ML, MR));
    const scalar fp = min(fp1Left, fp1Right);
    const scalar zeta = 1.0 - (1.0 - zeta1) * fp;

    const vector ULeft_star  = 0.5*(ULeft + URight_wall) + 0.5*zeta*(ULeft - URight_wall);
    const vector URight_star = 0.5*(ULeft + URight_wall) + 0.5*zeta*(URight_wall - ULeft);

    // --------------------------------------------------
    // 4. Conservative variables (reconstructed)
    // --------------------------------------------------
    const vector rhoULeft  = rhoLeft  * ULeft_star;
    const vector rhoURight = rhoRight * URight_star;

    const scalar rhoELeft =
        rhoLeft * (CvLeft * TLeft + 0.5 * magSqr(ULeft_star));

    const scalar rhoERight =
        rhoRight * (CvRight * TRight + 0.5 * magSqr(URight_star));

    const scalar qLeft  = ULeft_star  & n;
    const scalar qRight = URight_star & n;

    // --------------------------------------------------
    // 5. Roe averages (for wave speeds)
    // --------------------------------------------------
    const scalar rL = Foam::sqrt(max(0.0, rhoLeft));
    const scalar rR = Foam::sqrt(max(0.0, rhoRight));

    const scalar wL = rL / stabilise(rL + rR, VSMALL);
    const scalar wR = 1.0 - wL;

    const vector UTilde = wL * ULeft_star + wR * URight_star;
    const scalar unTilde = UTilde & n;

    const scalar HLeft  = (rhoELeft  + pLeft)  / rhoLeft;
    const scalar HRight = (rhoERight + pRight) / rhoRight;

    const scalar HTilde = wL * HLeft + wR * HRight;
    const scalar kappaTilde = wL * kappaLeft + wR * kappaRight;

    const scalar aTilde = Foam::sqrt
    (
        max(0.0, (kappaTilde - 1.0) * (HTilde - 0.5 * magSqr(UTilde)))
    );

    // --------------------------------------------------
    // 6. HLLE wave speeds (Einfeldt)
    // --------------------------------------------------
    const scalar SLeft  = min(qLeft - aLeft,  unTilde - aTilde);
    const scalar SRight = max(qRight + aRight, unTilde + aTilde);

    // --------------------------------------------------
    // 7. HLLE flux (NO star region!)
    // --------------------------------------------------
    const scalar invDen = 1.0 / stabilise(SRight - SLeft, VSMALL);

    // Left fluxes
    const scalar FL_rho = rhoLeft * qLeft;
    const vector FL_rhoU = rhoULeft * qLeft + pLeft * n;
    const scalar FL_rhoE = (rhoELeft + pLeft) * qLeft;

    // Right fluxes
    const scalar FR_rho = rhoRight * qRight;
    const vector FR_rhoU = rhoURight * qRight + pRight * n;
    const scalar FR_rhoE = (rhoERight + pRight) * qRight;

    // HLLE formula
    rhoFlux =
    (
        SRight * FL_rho
      - SLeft  * FR_rho
      + SRight * SLeft * (rhoRight - rhoLeft)
    ) * invDen * magSf;

    rhoUFlux =
    (
        SRight * FL_rhoU
      - SLeft  * FR_rhoU
      + SRight * SLeft * (rhoURight - rhoULeft)
    ) * invDen * magSf;

    rhoEFlux =
    (
        SRight * FL_rhoE
      - SLeft  * FR_rhoE
      + SRight * SLeft * (rhoERight - rhoELeft)
    ) * invDen * magSf;
}

// ************************************************************************* //
