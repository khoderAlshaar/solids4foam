
if(mesh.solutionDict().found("multiStage"))
{
    // Read multi-stage coeffs dict
    dictionary multiStage = mesh.solutionDict().subDict("multiStage");

    // Read number of subCycles
//     if (multiStage.found("numberSubCycles"))
//     {
//          multiStage.lookup("numberSubCycles") >> numberSubCycles;
//     }
//     if (multiStage.found("tolerance"))
//     {
//          multiStage.lookup("tolerance") >> tolerance;
//     }    
//     if (multiStage.found("relTol"))
//     {
//          multiStage.lookup("relTol") >> relTol;
//     }

}

