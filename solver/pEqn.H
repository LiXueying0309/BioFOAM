    volScalarField rAU("rAU", 1.0/UEqn.A()); // rAU:在速度的最后一个解中，从矩阵中提取对角项并存储倒数
//由于对流的非线性，矩阵系数是U的函数
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU)); 
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    volVectorField rhofHbyA ("rhofHbyA",rhof*HbyA);
//    volScalarField av ("av",mag(fvc::grad(eps)));
  //  volScalarField mdot ("mdot",R*av);
//UEqn.H是UEqn减去非对角项和U的乘积
//Jocobi 迭代
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rhof*rAU)*fvc::ddtCorr(U, phi)
    );
//将上述体积矢量场转化为面心标量场（参考有限体积法）
//fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)的作用是保证速度通量的全局守恒，以确保压力方程有解
    adjustPhi(phiHbyA, U, p);//调整边界值，保证速度边界条件守恒

    
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAUf);//更正压力边界条件，保证通量守恒
    
    while (piso.correctNonOrthogonal())// 网格非正交压力修正循环 
//如果在fvSolution字典文件中设置为0,就只求解控制方程一次,如果设置为n,则求解控制方程n-1次         
// Non-orthogonal pressure corrector loop
    {
        fvScalarMatrix pEqn// 构造压力方程进行求解
        (
           fvm::laplacian(rAUf,p) == fvc::div(phiHbyA)
               
        );

        pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));
//在不可压缩流动中，只有相对压力才起作用。除非存在边界压力
//一个网格的压力可以任意设定，以产生独特的压力解决方案

       pEqn.solve();

      if (piso.finalNonOrthogonalIter())//在最后一次非正交校正中，使用最新压力校正通量
        {
            phi = phiHbyA - pEqn.flux();

            p.relax();

           
	    //#include "correctUBc.H"
        }
   }
    #include "continuityErrs.H"//计算连续性方程误差

U = HbyA - rAU*fvc::grad(p);//校正速度，满足边界条件(主要针对第二类边界条件）
            U.correctBoundaryConditions();

    if (p.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
           pRefValue - getRefCellValue(p, pRefCell)
        );
    }

  

