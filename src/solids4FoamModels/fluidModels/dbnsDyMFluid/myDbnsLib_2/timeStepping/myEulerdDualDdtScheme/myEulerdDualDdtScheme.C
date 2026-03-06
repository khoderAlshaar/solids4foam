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

Class
    Foam::fv::backwardLocalDdtScheme

Author
    Oliver Borm
    Aleksandar Jemcov

\*---------------------------------------------------------------------------*/

#include "myEulerdDualDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar myEulerdDualDdtScheme<Type>::deltaT_() const
{
    return mesh().time().deltaT().value();
}



template<class Type>
scalar myEulerdDualDdtScheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0().value();
}


template<class Type>
template<class GeoField>
scalar myEulerdDualDdtScheme<Type>::deltaT0_(const GeoField& vf) const
{
    // Bug fix, Zeljko Tukovic: solver with outer iterations over a time-step
    // HJ, 12/Feb/2010
    if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
myEulerdDualDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt().internalField() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().V0() - coefft00*mesh().V00())/mesh().V()
        );

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
myEulerdDualDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*vf.internalField() -
                    (
                        coefft0*vf.oldTime().internalField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().internalField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
myEulerdDualDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.internalField() -
                    (
                        coefft0*vf.oldTime().internalField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().internalField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
myEulerdDualDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*rho.internalField()*vf.internalField() -
                    (
                        coefft0*rho.oldTime().internalField()
                       *vf.oldTime().internalField()*mesh().V0()
                      - coefft00*rho.oldTime().oldTime().internalField()
                       *vf.oldTime().oldTime().internalField()*mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*rho.boundaryField()*vf.boundaryField() -
                    (
                        coefft0*rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()
                      - coefft00*rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*rho*vf
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
myEulerdDualDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{

    // Info << "Hello from fvmDdt 1 "<<endl; 
    const objectRegistry& registry = this->mesh();


    const volScalarField& deltaTau =
        registry.lookupObject<volScalarField>(deltaTauName_);


    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();
    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);
    scalarField rDeltaTau = scalar(1.0)/(deltaTau.internalField());

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);              // = 3/2
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0)); // = 1/2
    scalar coefft0  = coefft + coefft00;                          // = 4/2
   
    // scalarField newTerm =  scalar(1.0)+ (scalar(2.0)*coefft *deltaTau.internalField() * rDeltaT);

    fvm.diag() = (rDeltaTau)*mesh().V();

    if (mesh().moving())
    {
        fvm.source() =
        rDeltaTau*vf.internalField()*mesh().V() - rDeltaT *
        (
              vf.internalField()*mesh().V()
            - vf.oldTime().internalField()*mesh().V0()
        );
    }

    else
    {
        fvm.source() =mesh().V() *
        (
            rDeltaTau*vf.internalField() - rDeltaT *
            (
                  vf.internalField()
                - vf.oldTime().internalField()
            )
        );
    }


    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
myEulerdDualDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const objectRegistry& registry = this->mesh();


    const volScalarField& deltaTau =
        registry.lookupObject<volScalarField>(deltaTauName_);

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);
    scalarField rDeltaTau = scalar(1.0)/(deltaTau.internalField());

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (rDeltaTau)*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rho.value()*
        (
            rDeltaTau*vf.internalField()*mesh().V() - rDeltaT*
            (
                  coefft*vf.internalField()*mesh().V()
                - coefft0*vf.oldTime().internalField()*mesh().V0()
            )
        );
    }

    else
    {
        fvm.source() =mesh().V()*rho.value()*
        (
            rDeltaTau*vf.internalField() - rDeltaT *
            (
                vf.internalField()
            - vf.oldTime().internalField()
            )
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
myEulerdDualDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
        Info << "Hello from fvmDdt 3 "<<endl; 

    const objectRegistry& registry = this->mesh();

    const volScalarField& deltaTau =
        registry.lookupObject<volScalarField>(deltaTauName_);

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);
    scalarField rDeltaTau = scalar(1.0)/(deltaTau.internalField());

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (rDeltaTau+coefft*rDeltaT)*rho.internalField()*mesh().V();

   fvm.diag() = (rDeltaTau)*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = 
              rDeltaTau*vf.internalField()*mesh().V()*rho.internalField()
            - rDeltaT*
            (
                  coefft*vf.internalField()*rho.internalField()*
                  mesh().V()
                - coefft0*vf.oldTime().internalField()*rho.oldTime().internalField()*
                  mesh().V0()
            );
    }

    else 
    {
        fvm.source() =mesh().V()*
        (
            rDeltaTau*vf.internalField()*rho.internalField() - rDeltaT *
            (
                coefft*vf.internalField()*rho.internalField()
            - coefft0*vf.oldTime().internalField()*rho.oldTime().internalField()
            )
        );
    }

    return tfvm;
}


template<class Type>
tmp<typename myEulerdDualDdtScheme<Type>::fluxFieldType>
myEulerdDualDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            rDeltaT*this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
           *(
                fvc::interpolate(rA)
               *(
                   coefft0*phi.oldTime()
                 - coefft00*phi.oldTime().oldTime()
                )
              - (
                    fvc::interpolate
                    (
                        rA*
                        (
                            coefft0*U.oldTime()
                          - coefft00*U.oldTime().oldTime()
                        )
                    ) & mesh().Sf()
                )
            )
        )
    );
}


template<class Type>
tmp<typename myEulerdDualDdtScheme<Type>::fluxFieldType>
myEulerdDualDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ','
      + rho.name() + ','
      + U.name() + ','
      + phiAbs.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT*this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime())
               *(
                    coefft0*fvc::interpolate(rA*rho.oldTime())
                   *phiAbs.oldTime()
                  - coefft00*fvc::interpolate(rA*rho.oldTime().oldTime())
                   *phiAbs.oldTime().oldTime()
                  - (
                        fvc::interpolate
                        (
                            rA*
                            (
                                coefft0*rho.oldTime()*U.oldTime()
                              - coefft00*rho.oldTime().oldTime()
                               *U.oldTime().oldTime()
                            )
                        ) & mesh().Sf()
                    )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                )
               *(
                    fvc::interpolate(rA*rho.oldTime())
                   *(
                       coefft0*phiAbs.oldTime()
                      /fvc::interpolate(rho.oldTime())
                     - coefft00*phiAbs.oldTime().oldTime()
                      /fvc::interpolate(rho.oldTime().oldTime())
                    )
                  - (
                        fvc::interpolate
                        (
                            rA*rho.oldTime()*
                            (
                                coefft0*U.oldTime()
                              - coefft00*U.oldTime().oldTime()
                            )
                        ) & mesh().Sf()
                    )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phiAbs.oldTime())
               *(
                    fvc::interpolate(rA)
                   *(
                       coefft0*phiAbs.oldTime()
                     - coefft00*phiAbs.oldTime().oldTime()
                    )
                  - (
                        fvc::interpolate
                        (
                            rA*
                            (
                                coefft0*U.oldTime()
                              - coefft00*U.oldTime().oldTime()
                            )
                        ) & mesh().Sf()
                    )
                )
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "myEulerdDualDdtScheme<Type>::fvcDdtPhiCorr"
        )   << "dimensions of phiAbs are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename myEulerdDualDdtScheme<Type>::fluxFieldType>
myEulerdDualDdtScheme<Type>::fvcDdtConsistentPhiCorr
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
)
{
    notImplemented
    (
        type()
      + "::fvcDdtConsistentPhiCorr"
      + "\n("
      + "\n    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,"
      + "\n    const GeometricField<Type, fvPatchField, volMesh>& U"
      + "\n    const surfaceScalarField rAUf"
      + "\n)"
    );

    // Dummy return
    return tmp<fluxFieldType>(nullptr);
}


template<class Type>
tmp<surfaceScalarField> myEulerdDualDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
        // Info << "Hello from myEulerdDualDdtScheme<Type>::meshPhi"<<endl;

    return mesh().phi();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
