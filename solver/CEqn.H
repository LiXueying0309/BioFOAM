//volScalarField av ("av",mag(fvc::grad(eps)));
//volScalarField k1 ("k1", Mmax*(mum/Yxs+ms));
//volScalarField k3 ("k3", Yxs/Mmax);
//volScalarField k4 ("k4", Mmax*ms);

 {   fvScalarMatrix CEqn
    (
	//Inertial Terms
        fvm::ddt(eps, CbyC0)
       + fvm::div(phi,CbyC0) 

	//Viscous Terms
      == fvm::laplacian(eps*D,CbyC0)
  
      - eps*Mmax/C0*(mum/Yxs+ms)*CbyC0*MbbyMmax/(ks/C0+CbyC0)

    );

    CEqn.solve();
}

    solve(UEqn==-fvc::grad(p));
//momentumPredictor: for low Reynolds number, off;
  if (piso.momentumPredictor()) 
    {
        solve
        (
            UEqn
            ==
            fvc::reconstruct
            (
               
                  - fvc::snGrad(p) * mesh.magSf()
		
            )
           
        );
    }
