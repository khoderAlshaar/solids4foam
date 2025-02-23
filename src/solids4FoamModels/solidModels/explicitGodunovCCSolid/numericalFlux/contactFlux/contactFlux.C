/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "contactFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace numericalFluxs
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(contactFlux, 0);
    addToRunTimeSelectionTable(numericalFlux, contactFlux, state);






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactFlux::contactFlux
(
    Time& runTime,
    const dynamicFvMesh& mesh,
    const volVectorField& lm,
    const volTensorField& F,
    const volTensorField& P,
    solidMaterialModel& model,
    operations& op,
    mechanics& mech
)
:
    numericalFlux(runTime),
    lm_(lm),
    F_(F),
    P_(P),
    model_(model),
    op_(op),
    mech_(mech),
    lmFlux_
    (
        IOobject
        (
            "lmFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(P_) & mesh.Sf())
    ),
    FFlux_
    (
        IOobject
        (
            "FFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       ((1/model_.density())*linearInterpolate(lm) * mesh.Sf())
    )

{

    Info << "Hello from contactFlux constructor" << endl;
   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
contactFlux::~contactFlux()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void contactFlux::computeFlux()
{
    Info << "Compute flux using Contact flux" << endl;
    // Info << lm_ <<endl;
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numericalFluxs

} // End namespace Foam

// ************************************************************************* //
