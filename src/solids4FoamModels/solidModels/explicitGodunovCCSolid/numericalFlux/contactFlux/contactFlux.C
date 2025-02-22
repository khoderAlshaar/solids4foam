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
    Time& runTime
)
:
    numericalFlux(runTime)   
    // numericalFlux()   
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
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numericalFluxs

} // End namespace Foam

// ************************************************************************* //
