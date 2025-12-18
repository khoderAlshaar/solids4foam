#include "externalStiffenedGasThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(externalStiffenedGasThermo, 0);

    addToRunTimeSelectionTable
    (
        basicPsiThermo,
        externalStiffenedGasThermo,
        fvMesh
    );
}


// ===== Internal calculations =====
void Foam::externalStiffenedGasThermo::calculate()
{

    const scalarField& hCells = h_.internalField();
    const scalarField& pCells = this->p_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();

    scalarField& rhoCells = rho_.internalField();


    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();


    forAll(TCells, celli)
    {

        TCells[celli] = TH(hCells[celli], TCells[celli]);

        rhoCells[celli] = rho(pCells[celli], TCells[celli]); 

        psiCells[celli] = psi(pCells[celli], TCells[celli]);

        muCells[celli] = mu(TCells[celli]);
        alphaCells[celli] = alpha(TCells[celli]);

    } 

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];

        fvPatchScalarField& ph = h_.boundaryField()[patchi];

        fvPatchScalarField& prho = rho_.boundaryField()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];

        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];


        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                ph[facei] = H(pT[facei]); //! if T is fixd value update h from T is this right? 
                                         //! but h is already updated in updateFields.H from conservative variables 
                // pT[facei] = TH(ph[facei], pT[facei]);

                prho[facei] = rho(pp[facei], pT[facei]);

                ppsi[facei] = psi(pp[facei], pT[facei]);


                pmu[facei] = mu(pT[facei]);

                palpha[facei] = alpha(pT[facei]);


            }
        }
        else
        {
            forAll(pT, facei)
            {
                pT[facei] = TH(ph[facei], pT[facei]);

                prho[facei] = rho(pp[facei], pT[facei]);

                ppsi[facei] = psi(pp[facei], pT[facei]);

                pmu[facei] = mu(pT[facei]);
                
                palpha[facei] = alpha(pT[facei]);
               
            }
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalStiffenedGasThermo::externalStiffenedGasThermo
(
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    basicPsiThermo(mesh, obj),

    mesh_(mesh),

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            obj,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->hBoundaryTypes()
    ),
    
    rho_
    (
        IOobject
        (
            "rhoSGThermo",
            mesh.time().timeName(),
            obj,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    gamma_(0.0),
    pInf_(0.0),
    Cp_(0.0),
    Cv_(0.0),
    R_(0.0),
    muRef_(0),
    Pr_(0),
    q_(0)
{

    Info <<"Hello from externalStiffenedGasThermo constructor"<<endl;
    read();
    Cp_ = Cv_*gamma_;
    R_ = (gamma_- 1)*Cv_;
    // Cv_ = Cp_/(gamma_ - 1);
    // R_  = Cp_ - Cv_;

    scalarField& hCells = h_.internalField();
    const scalarField& TCells = this->T_.internalField();
    
    forAll(hCells, celli)
    {
        hCells[celli] = H(TCells[celli]);
    }

    forAll(h_.boundaryField(), patchi)
    {
        h_.boundaryField()[patchi] =
            h(this->T_.boundaryField()[patchi], patchi);

    }

    hBoundaryCorrection(h_);
    

    calculate();
}


// ===== Destructor =====

Foam::externalStiffenedGasThermo::~externalStiffenedGasThermo()
{}



// ===== Scalar Evaluations =====

Foam::tmp<Foam::scalarField>
Foam::externalStiffenedGasThermo::h
(
    const scalarField& T,
    const labelList&
) const
{
    return gamma_*Cv_*T + q_;
}

Foam::tmp<Foam::scalarField>
Foam::externalStiffenedGasThermo::h
(
    const scalarField& T,
    const label
) const
{
    return gamma_*Cv_*T + q_;
}


// ===== Cp =====

Foam::tmp<Foam::scalarField>
Foam::externalStiffenedGasThermo::Cp
(
    const scalarField& T,
    const labelList&
) const
{
    return tmp<scalarField>(new scalarField(T.size(), Cp_));
}

Foam::tmp<Foam::scalarField>
Foam::externalStiffenedGasThermo::Cp
(
    const scalarField& T,
    const label
) const
{
    return tmp<scalarField>(new scalarField(T.size(), Cp_));
}

Foam::tmp<Foam::volScalarField>
Foam::externalStiffenedGasThermo::Cp() const
{

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Cp", dimensionSet(0, 2, -2, -1, 0), Cp_)
        )
    );

    return tCp;
}


// ===== Cv =====

Foam::tmp<Foam::scalarField>
Foam::externalStiffenedGasThermo::Cv
(
    const scalarField& T,
    const labelList&
) const
{
    return tmp<scalarField>(new scalarField(T.size(), Cv_));
}

Foam::tmp<Foam::scalarField>
Foam::externalStiffenedGasThermo::Cv
(
    const scalarField& T,
    const label
) const
{
    return tmp<scalarField>(new scalarField(T.size(), Cv_));
}

Foam::tmp<Foam::volScalarField>
Foam::externalStiffenedGasThermo::Cv() const
{
    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Cv", dimensionSet(0, 2, -2, -1, 0), Cv_)
        )
    );

    return tCv;

}



// ===== Correction =====

void Foam::externalStiffenedGasThermo::correct()
{
    calculate();
}


// ===== Dictionary =====

bool Foam::externalStiffenedGasThermo::read()
{
    if (basicPsiThermo::read())
    {
        // Enclose the creation of the thermophysicalProperties to ensure it is
        // deleted before the turbulenceModel is created otherwise the dictionary
        // is entered in the database twice
        {
            IOdictionary thermoDict
            (
                IOobject
                (
                    "thermophysicalProperties",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            );

            const dictionary& d = thermoDict.subDict("stiffenedGasCoeffs");
            d.lookup("gamma") >> gamma_;
            d.lookup("pInf") >> pInf_;
            // d.lookup("Cp") >> Cp_;
            d.lookup("Cv") >> Cv_;
            d.lookup("q") >> q_;
            d.lookup("mu") >> muRef_;
            d.lookup("Pr") >> Pr_;
        }




        return true;
    }




    return false;
}
