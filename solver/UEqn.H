    fvVectorMatrix UEqn//构造体积矢量场矩阵
    (
	//Inertial Terms
        (1./eps)*(fvm::ddt(rhof, U) + fvm::div(rhoPhiByEps, U) )

	//Viscous Terms
      - (1./eps)*(fvm::laplacian(muf,U))
      

	//Implicit Darcy-Brinkman Drag Term
     + fvm::Sp(Drag,U) 
 

    );

    UEqn.relax();
 // solve(UEqn==-fvc::grad(p));
//如果进行速度预测（on），则求解完整的动量方程得到预测速度
//如果不进行速度预测(off），则预测速度直接取当前已知时间步的速度
//momentumPredictor：动量预测求解开关，对多相流以及低雷诺数一般设置为off；
 if (piso.momentumPredictor()) 
    {
        solve
        (
            UEqn
           == -fvc::grad(p)
        );
    }

