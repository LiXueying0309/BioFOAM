{
    fvScalarMatrix epsEqn
    (
       fvm::ddt(eps)  ==  -(1./rhoMb)*fvc::ddt(MbbyMmax)
    );

	epsEqn.solve();
	eps.correctBoundaryConditions();
	eps.min(0.00000001);
	eps.max(1);
}

forAll(eps,celli) 
{
	if(eps[celli]>0.5) // 0.5 is a sample value for the adjustable pore volume fraction within a valid range. 
	{	
		Solid[celli]=0; 
	}
	else
	{
		Solid[celli]=1; 
	}
}

