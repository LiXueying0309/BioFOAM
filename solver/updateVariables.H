
volScalarField k0("k0", K*(eps*eps*eps)/(max((1-eps)*(1-eps),SMALL)) ); //Uncomment for Kozeny-Carman
//volScalarField k0("k0", K); 
surfaceScalarField k0f("k0f",fvc::interpolate(k0,"k0"));


//Mobilities
volScalarField M ("M",k0/muf); 
//volScalarField L1 ("L",rho1*k0*kr1/mu1);	


//Drag Coefficient Calculation
//surfaceScalarField Dragf ("Dragf",1/Mf);
volScalarField Drag ("Drag", 1/M);

forAll(Drag,celli)
{
	if (Solid[celli] == 0)
		{Drag[celli]= 0;}

}
