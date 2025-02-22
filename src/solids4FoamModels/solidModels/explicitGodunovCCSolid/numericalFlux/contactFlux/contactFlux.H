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

Class
    contactFlux

Description
    Virtual base class for solid mechanics models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    contactFlux.C
    contactFluxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef contactFlux_H
#define contactFlux_H

#include "numericalFlux.H"
// #include "IOdictionary.H"
// #include "scalarIOField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

namespace numericalFluxs
{
/*---------------------------------------------------------------------------* \
                           Class contactFlux Declaration
\*---------------------------------------------------------------------------*/

class contactFlux
:
    public numericalFlux
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        contactFlux(const contactFlux&);

        //- Disallow default bitwise assignment
        void operator=(const contactFlux&);


protected:

    // Protected member functions


public:

    //- Runtime type information
    TypeName("contact");


    // Constructors

        //- Construct from components
        contactFlux
        (
            Time& runTime
        );

    // Destructor

        virtual ~contactFlux();


    // Member Functions

        // Access

        // Edit

            //- Evolve the solid model
            virtual void computeFlux();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numericalFluxs


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
