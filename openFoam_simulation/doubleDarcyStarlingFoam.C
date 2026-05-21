#include "fvCFD.H"
#include "pimpleControl.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createTimeControls.H"

    Info<< "Reading fields: p1, p2, K1, K2, Chi2, C1, C2\n" << endl;

    volScalarField p1(IOobject("p1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField p2(IOobject("p2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);

    volScalarField Chi2(IOobject("Chi2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);

    volScalarField K1(IOobject("K1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField K2(IOobject("K2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);

    volScalarField C1(IOobject("C1", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField C2(IOobject("C2", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);

    IOdictionary props
    (
        IOobject("doubleDarcyProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
    );

    dimensionedScalar mu("mu", props);
    dimensionedScalar Lp("Lp", props);
    dimensionedScalar beta("beta", props);
    dimensionedScalar D1("D1", props);
    dimensionedScalar D2("D2", props);
    dimensionedScalar P_DL("P_DL", props);
    dimensionedScalar S_mLV("S_mLV", props);
    dimensionedScalar V_dura("V_dura", props);
    dimensionedScalar phi_D("phi_D", props);
    dimensionedScalar phi_L("phi_L", props);
    dimensionedScalar R("R", props);
    dimensionedScalar Temperature("Temperature", props);
    dimensionedScalar sigma("sigma", props);



    // --- CHI2 SMOOTHING (HELMHOLTZ FILTER) ---
    Info<< "\nApplying Helmholtz filter to smooth Chi2 interface...\n" << endl;
.
    dimensionedScalar lambda2("lambda2", dimArea, 5e-3); // Start with 1e-5, increase if you need more blur

    // Store the original sharp 0/1 field as the source term
    volScalarField Chi2_sharp
    (
        IOobject("Chi2_sharp", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        Chi2
    );

    // Solve the implicit spatial filter
    fvScalarMatrix Chi2FilterEqn
    (
        fvm::Sp(1.0, Chi2)
      - fvm::laplacian(lambda2, Chi2)
     == Chi2_sharp
    );

    Chi2FilterEqn.solve();

    // Clamp the field to prevent any numerical overshoots beyond physical bounds
    Chi2 = max(Chi2, dimensionedScalar("zero", dimless, 0.0));
    Chi2 = min(Chi2, dimensionedScalar("one", dimless, 1.0));

    Chi2.write();

    Info<< "Chi2 smoothing complete.\n" << endl;
    // ----------------------------------------


    // --- STEADY-STATE PRESSURE SOLVE ---
    Info<< "\nSolving steady-state pressure fields (Homotopy stabilization)...\n" << endl;

    int maxOuterIter = 2000;
    int nNonOrthCorr = 3;
    scalar tolerance = 1e-4;

    // --- Relax (Homotopy) ---
    scalar dpScaleStart = 100.0;
    scalar dpScaleMin = 0.01;
    scalar iterDecay = 30.0;

    volScalarField smoothPos
	(
    		IOobject("smoothPos", runTime.timeName(), mesh,IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar(dimless, 0.0)
);

    for (int iter = 0; iter < maxOuterIter; ++iter)
    {


	scalar currentDpScale = dpScaleMin + (dpScaleStart - dpScaleMin) * Foam::exp(-scalar(iter) / iterDecay);

    	smoothPos = 0.5 * (1.0 + tanh((p1 - p2) / dimensionedScalar("dp", dimPressure, currentDpScale)    ));

        scalar p1Res = 0.0;
        scalar p2Res = 0.0;

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; ++nonOrth)
        {
            fvScalarMatrix p1Eqn
            (
                -fvm::laplacian(K1/mu, p1)
                == smoothPos*Chi2*fvm::Sp(-Lp*S_mLV/V_dura, p1)
                 + smoothPos*Chi2*Lp*(S_mLV/V_dura)*p2
            );

            SolverPerformance<scalar> solverPerf_p1 = p1Eqn.solve();

            fvScalarMatrix p2Eqn
            (
                -fvm::laplacian(K2/mu, p2)
                == smoothPos*fvm::Sp(-Chi2*Lp*S_mLV/V_dura, p2)
                 + smoothPos*Chi2*Lp*(S_mLV/V_dura)*p1
            );

            SolverPerformance<scalar> solverPerf_p2 = p2Eqn.solve();

            if (nonOrth == 0)
            {
                p1Res = solverPerf_p1.initialResidual();
                p2Res = solverPerf_p2.initialResidual();
            }
        }

        // Print
        Info<< "Iter: " << iter
            << " | dpScale: " << currentDpScale
            << " | p1 Res: " << p1Res
            << " | p2 Res: " << p2Res << endl;

        if (p1Res < tolerance && p2Res < tolerance && currentDpScale < (dpScaleMin * 1.1))
        {
            Info<< "\nSteady-state pressure converged dynamically in " << iter << " iterations." << endl;
            break;
        }

        if (iter == maxOuterIter - 1)
        {
            WarningInFunction
                << "Steady-state pressure solve hit max iterations (" << maxOuterIter
                << ") without fully converging." << endl;
        }
    }

    // Steady Flux
    surfaceScalarField phi1 = linearInterpolate(-(K1/mu)*fvc::grad(p1)) & mesh.Sf();
    surfaceScalarField phi2 = linearInterpolate(-(K2/mu)*fvc::grad(p2)) & mesh.Sf();


    // Pressure heaviside
    volScalarField posCoupling = Foam::pos(p1 - p2);

    surfaceScalarField phi1_norm = phi1 / phi_D;
    surfaceScalarField phi2_norm = phi2 / phi_L;

    volScalarField pi1
    (
        IOobject("pi1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        R * Temperature * C1
    );
    volScalarField pi2
    (
        IOobject("pi2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        R * Temperature * C2
    );

    Info<< "Steady-state pressure solve complete.\n" << endl;
    // -----------------------------------

    Info<< "\nStarting Time Loop for Scalar Transport\n" << endl;
    volScalarField couplingC1toC2
    (
        IOobject("couplingC1toC2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        posCoupling * Chi2 * P_DL
    );

    volScalarField couplingC2toC1
    (
        IOobject("couplingC2toC1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        couplingC1toC2 * (phi_D / phi_L)
    );

    volScalarField starlingFlux
    (
        IOobject("starlingFlux", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        posCoupling * Chi2 * Lp * (S_mLV/V_dura) * (p2 - p1 - sigma*(pi1 - pi2)) * 0.5
    );
    while (runTime.run())
    {
        #include "readTimeControls.H"

        // Courant calculation for AdjustTime (If Needed)
        scalar CoNum1 = 0.0;
        scalar CoNum2 = 0.0;
        scalar trueCoNum1 = 0.0;
        scalar trueCoNum2 = 0.0;
        scalar meanCoNum = 0.0;

        if (mesh.nInternalFaces() > 0)
        {
            scalarField sumPhi1(fvc::surfaceSum(mag(phi1))().internalField());
            scalarField sumPhi2(fvc::surfaceSum(mag(phi2))().internalField());

            CoNum1 = 0.5 * gMax(sumPhi1 / mesh.V().field()) * runTime.deltaTValue();
            CoNum2 = 0.5 * gMax(sumPhi2 / mesh.V().field()) * runTime.deltaTValue();

            trueCoNum1 = CoNum1 / phi_D.value();
            trueCoNum2 = CoNum2 / phi_L.value();

            scalar meanTrueCo1 = (0.5 * (gSum(sumPhi1) / gSum(mesh.V().field())) * runTime.deltaTValue()) / phi_D.value();
            scalar meanTrueCo2 = (0.5 * (gSum(sumPhi2) / gSum(mesh.V().field())) * runTime.deltaTValue()) / phi_L.value();
            meanCoNum = max(meanTrueCo1, meanTrueCo2);
        }

        scalar CoNum = max(trueCoNum1, trueCoNum2);

        Info<< "Bulk Courant max: p1=" << CoNum1 << ", p2=" << CoNum2 << nl
            << "True Scalar Courant max: C1=" << trueCoNum1 << ", C2=" << trueCoNum2 << nl
            << "Overall Max Courant Number: " << CoNum << " Mean: " << meanCoNum << endl;

        #include "setDeltaT.H"

        runTime++;

        // C1 and C2 BCs update in the time loop
        C1.correctBoundaryConditions();
        C2.correctBoundaryConditions();

        pi1 = R * Temperature * C1;
        pi2 = R * Temperature * C2;

        couplingC1toC2 = posCoupling * Chi2 * P_DL;
        couplingC2toC1 = couplingC1toC2 * (phi_D / phi_L);
        starlingFlux   = posCoupling * Chi2 * Lp * (S_mLV/V_dura) * (p2 - p1 - sigma*(pi1 - pi2)) * (1-sigma)*0.5;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix C1Eqn
                (
                    fvm::ddt(C1)
                  + fvm::div(phi1/phi_D, C1)
                  - fvm::laplacian(D1,C1)
                 ==couplingC2toC1 * C2
                  + fvm::Sp(-couplingC1toC2, C1)
                  + fvm::Sp(starlingFlux/phi_D, C1)

                );
                C1Eqn.relax();
                C1Eqn.solve();
                //C1.clamp_min(0.0);

                fvScalarMatrix C2Eqn
                (
                    fvm::ddt(C2)
                  + fvm::div(phi2/phi_L, C2)
                  - fvm::laplacian(Chi2*D2,C2)
                 ==couplingC1toC2 * C1
                  + fvm::Sp(-couplingC2toC1, C2)
                  - (starlingFlux/phi_D) * C1
                );
                C2Eqn.relax();
                C2Eqn.solve();

            } // EOW correctNonOrthogonal
        } // EOW pimple

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime()
            << " s  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

    } // EOW time

    Info<< "End\n" << endl;
    return 0;
}
