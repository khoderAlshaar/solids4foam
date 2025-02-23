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

// Strain measures
H_ = det(F_)*op_.invT(F_);
J_ = det(F_);

// Wave speeds
mech_.correct(Up_, Us_, F_);
Up_time_ = model_.Up()/mech_.stretch();

// Constitutive model
model_.correct();
P_ = model_.piola();

flux().computeFlux();
