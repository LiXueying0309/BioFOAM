
volScalarField k0("k0", K*(max((1-eps)*(1-eps),SMALL))/(eps*eps*eps)); //Uncomment for Kozeny-Carman
surfaceScalarField k0f("k0f",fvc::interpolate(k0,"k0"));

//Mobilities
volScalarField M ("M",k0/muf); 
volScalarField Drag ("Drag", 1/M);
forAll(Drag,celli)
{
	if (Solid[celli] == 0)
		{Drag[celli]= 0;}

}
