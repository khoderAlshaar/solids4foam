#ifndef externalStiffenedGasThermo_H
#define externalStiffenedGasThermo_H

#include "basicPsiThermo.H"
#include "volFields.H"
#include "IOdictionary.H"

namespace Foam
{

class externalStiffenedGasThermo
:
    public basicPsiThermo
{
    // Private data

        const fvMesh& mesh_;

        volScalarField h_;
        volScalarField rho_;

        scalar gamma_;
        scalar pInf_;
        scalar Cp_;
        scalar Cv_;
        scalar R_;
        scalar muRef_;   // scalar reference viscosity (used in alpha())
        scalar Pr_;
        scalar q_;


    // Private member functions

        void calculate();

        inline scalar T
        (
            scalar f,
            scalar T0
        ) const;

        // Disallow copy
        externalStiffenedGasThermo
        (
            const externalStiffenedGasThermo&
        );


public:

    TypeName("externalStiffenedGasThermo");


    //  REQUIRED CONSTRUCTOR
    externalStiffenedGasThermo
    (
        const fvMesh& mesh,
        const objectRegistry& obj
    );



    // Destructor
    virtual ~externalStiffenedGasThermo();


    // Member functions

        //- Enthalpy [J/kg]
        inline scalar H(const scalar T) const;

        // Energy->temperature  inversion functions
        //- Temperature from Enthalpy given an initial temperature T0
        inline scalar TH(const scalar H, const scalar T0) const;
       
        // - Thermal diffusivity for enthalpy [kg/ms]
        inline scalar alpha(const scalar T) const;
        // Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar) const;

                //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar h) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

     //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }


        // Fields derived from thermodynamic state variables
            //- Density [kg/m^3] 
            //! this one should be called only at creatFields in the solver for initiazation
            //! and to disallow the function in basicPsiThermo.H in the case of stiffend gas EOS 
            virtual tmp<volScalarField> rho() const
            {
                return rho_;
            }



            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for cell-set [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();


        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "externalStiffenedGasThermoI.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
