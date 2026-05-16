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

Class
    SGTNRoeFlux

Description
    Implementation of inviscid numerical fluxes based on Roe & Pike FDS scheme
    with Harten's entropy fix. This is a first order scheme

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    SGTNRoeFlux.H

\*---------------------------------------------------------------------------*/

#ifndef SGTNRoeFlux_H
#define SGTNRoeFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SGTNRoeFlux Declaration
\*---------------------------------------------------------------------------*/

class SGTNRoeFlux:
public dbnsFlux
{
    // privet 
    scalar pInf;
    scalar gamma;
    scalar Cv;
    scalar q;
public:
    TypeName("SGTNRoeFlux");

    // Constructors

    //- Construct null
    SGTNRoeFlux(const fvMesh& mesh, const dictionary&)
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
                d.lookup("pInf") >> pInf;
                d.lookup("gamma") >> gamma;
                d.lookup("Cv") >> Cv;
                d.lookup("q") >> q;
            }
            else
            {
                FatalErrorIn("SGTNRoeFlux")
                << "thermo Type Name = " << thermoTypeName << nl
                << " it shoudl be externalStiffenedGasThermo " << nl
                << exit(FatalError);
            }
        }
    
    
    // Destructor - default
    virtual ~SGTNRoeFlux() {};

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
