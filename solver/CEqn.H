//volScalarField av ("av",mag(fvc::grad(eps)));
//volScalarField k1 ("k1", Mmax*(mum/Yxs+ms));
//volScalarField k3 ("k3", Yxs/Mmax);
//volScalarField k4 ("k4", Mmax*ms);

 {   fvScalarMatrix CEqn//构造体积矢量场矩阵
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
//如果进行速度预测（on），则求解完整的动量方程得到预测速度
//如果不进行速度预测(off），则预测速度直接取当前已知时间步的速度
//momentumPredictor：动量预测求解开关，对多相流以及低雷诺数一般设置为off；
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

