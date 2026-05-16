/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
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

Application
    hllcSGALEFlux

Description
    HLLC flux class, computes numerical fluxes based on HLLC scheme. Uses
    primitive variables as input and gives back conservative numerical fluxes.
    
    Stiffened equation of state is used to compute the required variables from 
    given input.
    
    Toro, E.I. and Spruce, M. and Spears, W.: Restoration of the Contact
    Surface in the HLL-Riemann Solver, Shock Waves, 4:25-34, 1996
    
Author
    Ghanshyam Bharate, IIT Bombay

SourceFiles
    hllcSGALEFlux.H

\*---------------------------------------------------------------------------*/

#ifndef hllcSGALEFlux_H
#define hllcSGALEFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hllcFlux Declaration
\*---------------------------------------------------------------------------*/

class hllcSGALEFlux:
public dbnsFlux
{
    // privet member data
    scalar pinf;
    scalar gamma;
    scalar Cv;
    scalar q;

public:
    TypeName("hllcSGALEFlux");
    // Constructors

        //- Constructor null
        hllcSGALEFlux(const fvMesh& mesh, const dictionary&)
        {
            word thermoTypeName;
            
                IOdictionary thermoDict
                (
                    IOobject
                    (
                        "thermophysicalProperties",
                        mesh.time().constant(),
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                thermoDict.lookup("thermoType") >> thermoTypeName;
            
            if (thermoTypeName == "externalStiffenedGasThermo")
            {
                const dictionary& d = thermoDict.subDict("stiffenedGasCoeffs");
                d.lookup("pInf") >> pinf;
                d.lookup("gamma") >> gamma;
                d.lookup("Cv") >> Cv;
                d.lookup("q") >> q;
            }

            // if (thermoTypeName == "eRhoThermo<constTransportEos<SGEOS>>")
            // {
                // const dictionary& d = thermoDict.subDict("stiffenedGasCoeffs");
                // d.lookup("pInf") >> pInf;
                // d.lookup("gamma") >> gamma;
                // d.lookup("Cv") >> Cv;
                // d.lookup("q") >> q;
                 
                // Read thermodynamic properties from thermo dictionary
                // word firstToken;
                // thermoDict.lookup("fluidProperties") >> firstToken >> q
                // >> pinf >> gamma >> Cv;
            // }
            else
            {
                FatalErrorIn("hllcSGALEFlux")
                << "thermo Type Name = " << thermoTypeName << nl
                << " it shoudl be eRhoThermo " << nl
                << exit(FatalError);
            }
        }


    // Destructor - default
     virtual ~hllcSGALEFlux() {};


    // Member Functions

        //- Update fluxes based on current left and right states
      virtual void evaluateFlux
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
         ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
