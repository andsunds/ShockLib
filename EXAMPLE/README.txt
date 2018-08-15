**********************************************************************
************* IMPORTANT NOTE *****************************************
**********************************************************************
You must have an enviroment variabe
   SHOCKLIB
leading to the Shock library, e.g.
/home/andsunds/SVN/erc/andsunds/Shock-project/ShockLib/

This is done by executing:
   export SHOCKLIB=<path to ShockLib directory>
in the terminal, or adding the above line to your .bashrc file
(or equivalent).
**********************************************************************
**********************************************************************
**********************************************************************

This directory should contain a piece of (almost) self-contained code for
generating analytical approximations of out semi-analytical shocks and
importing them into Gkyl.

This is done through an example, found in the directory: EXAMPLE/
1. shock_creation_EX.m	- The first step is to run through the shock
   		       	  creation script with the desired parameters.
			  This produces a "psi-E-rho_<...>.tsv" file
			  with raw data.
2. shock_approx_EX.m	- Next, run the script for finding the shock
   		       	  approximation, tweaking the parameters until
			  saftisfied. This produces a "COEFS-<...>.tsv"
			  file with the coefficients.
3. g2-Shock-EX.lua	- Finally, you can run the Gkyl script, with the
   			  desired "COEFS-<...>.tsv" file as the input,
			  and the desired simulation paramters.

The idea/hope is that this is a relatively well packaged code package.
In particular, it should be possible to copy the whole EXAMPLE/ folder
to somewhere more convenient, everything should still work.


There will also be a more detailed instruction in the instrucion notes
of the whole Shock library.


