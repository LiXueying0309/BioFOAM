volTensorField gradU(fvc::grad(U));
volScalarField shearStress = muf * gradU.component(tensor::XY);
volScalarField phi_det = pos(shearStress - tau_max) * (1 - exp(-C_det * (shearStress - tau_max) / tau_max));
// The detachment rate might be a function of phi_det and the biomass density
volScalarField detachmentRate = phi_det * detachmentValue;
volScalarField attachmentRate = pos(U.component(0) - U_attach) * attachmentValue;

{
volScalarField conditionalAttachmentRate
        (
            IOobject
            (
                "conditionalAttachmentRate",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0) // Initialize with 0
        );

        forAll(mesh.C(), celli) // Loop over all cells
        {
            if (mesh.C()[celli].x() < 1.25 && mesh.C()[celli].x() > -1.25) // Check x-coordinate of cell center, convert mm to meters
            {
                conditionalAttachmentRate[celli] = attachmentRate[celli];
            }
        }

fvScalarMatrix MEqn
(
    fvm::ddt(eps, MbbyMmax) 
    ==
    eps*fvm::laplacian(Db, MbbyMmax)
    + fvm::div(phi, MbbyMmax)
    + eps*(mum/Yxs + ms)*Yxs*CbyC0*MbbyMmax/(ks/C0 + CbyC0)
    - eps*Yxs*ms*MbbyMmax
    + eps*conditionalAttachmentRate
    - eps*detachmentRate
);

MEqn.solve();

//#include "epsEqn.H"

    forAll(Solid,celli) 
    {
      if(Solid[celli]==1) 
      
                { //MbbyMmax[celli] = 1; 
                MbbyMmax[celli] =0;}
    }
  
   forAll(MbbyMmax,celli) 
    {
      if(MbbyMmax[celli]>0.3) 
      
                { MbbyMmax[celli] = 0.3; 
                eps[celli] =0.0001;}
      if(MbbyMmax[celli]<-0.0001) 
      
                { MbbyMmax[celli] = 0; 
                eps[celli] =1;}
    }


 /*  forAll(MbbyMmax,celli)
   {
       if(MbbyMmax[celli] >= 1)
        {Solid[celli] = 1;
         eps[celli] =0;
        }
        else
        {Solid[celli] = 0;
        eps[celli] = 1;

        }


        forAll(eps,celli) 
   {
    if(eps[celli]>0.1) 
         {
            if(MbbyMmax[celli]<1) 
      
                {MEqn.solve(); }
  }


  else
    {
    MbbyMmax[celli] = 0;
    }
*/



   }
 
 //  MbbyMmax.max(1);
//cout << "BEGIN" << endl;
//cout << "cellCells():" << cellCells() << endl;
//cout << "BEGIN" << endl;


//if ( mesh().lookupObject<volScalarField>("Mb") > 4 )
//{ eps=0.000001; }
/*
forAll(MbbyMmax, celli)
{
    if (MbbyMmax[celli] > 1)
    {
        eps[celli] = 0;
      //  Solid[celli]=1;
       U[celli].component(0) = 0;
       U[celli].y()=0;
       U[celli].component(2)=0;
    }
    forAll(eps,celli) 
    {
  	if(eps[celli]>0.99) 
	{	
	Solid[celli]=0; 
	}
	else
	{
	Solid[celli]=1; 
	}
    }

}

*/



