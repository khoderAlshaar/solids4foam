/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "externalStiffenedGasThermo.H"

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //


namespace Foam
{

inline Foam::scalar Foam::externalStiffenedGasThermo::T
(
    scalar f,
    scalar T0
) const
{
    // Info << "TEMPRETURE is evaluated in speciesThermoL.H SECOND" <<endl;
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*1e-4;
    int iter = 0;
    // Info<< "TH(),h= " << f <<endl;
    do
    {
        // Limit the temperature jump in a single corrector to TJump_
        // HJ, 12/Oct/2010

        Test = Tnew;
        Tnew = Test
            - Foam::min((H(Test) - f)/Cp_, 20);

        if (iter++ > 100)
        {
            // Improvements: graceful exit with recovery.  HJ, 11/Oct/2010
            InfoIn
            (
                "externalStiffenedGasThermo::T(scalar f, scalar T0, "
                "scalar (externalStiffenedGasThermo::*F)(const scalar) const, "
                "scalar (externalStiffenedGasThermo::*dFdT)(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded.  T0 = "
                << T0 << " F = " << H(Test)
                << " dFdT = " << Cp_
                << endl;

            // Use value where dFdT is calculated using T0.  HJ, 11/Oct/2010
            Tnew = f/Cp_;

            return Tnew;
        }

        // Info << "Tnew - Test = "<<Tnew - Test<<endl;

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}



inline Foam::scalar Foam::externalStiffenedGasThermo::H(const scalar T) const
{
    // Info << "ENTHALPY is evaluated in stiffenedGasSpecieThermo.H from H(const scalar T) " <<endl;
    return  gamma_*Cv_ * T + q_;
}


inline Foam::scalar Foam::externalStiffenedGasThermo::TH
(
    const scalar h,
    const scalar T0
) const
{
    // Info << "TEMPRETURE is evaluated in speciesThermoL.H FIRST" <<endl;
    // return T(h, T0);
    return (h - q_)/(gamma_*Cv_);
}



// ===== Transport =====
inline Foam::scalar  Foam::externalStiffenedGasThermo::mu(const scalar) const
{
    return muRef_;
}

inline  Foam::scalar  Foam::externalStiffenedGasThermo::alpha(const scalar T) const
{
    scalar Tstd = 298.15; // reference temperature

    scalar deltaT = T - Tstd;
    scalar CpBar =
        (deltaT*(H(T) - H(Tstd)) + Cp_)/(sqr(deltaT) + 1);
    scalar rPr = 1.0/Pr_;

    // use scalar reference viscosity (muRef_) rather than field mu_
        // return Cp_*muRef_*rPr / CpBar;
        return Cp_*muRef_*rPr;

    // return muRef_/(rho*Pr_);
}


// ===== Equation of State =====

inline  Foam::scalar externalStiffenedGasThermo::rho(scalar p, scalar T) const
{
    // Info << "Density is evaluated in stiffenedGasL.H"<<endl;
    return (p + pInf_)/((gamma_ - 1)*Cv_ *T) ;
}


//- Return compressibility rho/p [s^2/m^2]
inline Foam::scalar externalStiffenedGasThermo::psi(const scalar p, const scalar rho) const
{
    scalar c2 = gamma_*(p + pInf_)/rho;
    return 1.0/(rho*c2);
    // return 1.0/(R_*T);

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
