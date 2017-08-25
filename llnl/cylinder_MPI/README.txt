Cylinder to ParaDiSv2.5.1 Comments
----------------------------------------------------------------------------------

Wei Cai Sep 15 2011

This version is created by convert_thinfilm_to_cylinder








Cylinder from ParaDiS to ParaDiSv2.2.6 Comments
----------------------------------------------------------------------------------

Sylvie Nov 23 2009

- Removed reference to cylinder0 and home0
- Modified ParadisStep.c and ParaDiSProto.h
- Straighten out PrinStress routine
- Changed CYL_Main.c, CYL_Util.c, cylinder.c to mimic thinfilm code.
  tractions are computed wilt AllSegmentStress, 3x3 stresses, MPI_reduce
  on tractions and not on stress



To do list:
- Transform cart2cyl and cyl2cart from 6 stresses to 3x3 stresses
  It is a bit tricky

-Need a lot of checks!!!!!



--------------------------------------------------------------------------------

Wei.   Wed, Aug 5 2009

1. Initialize.c, copied from src/ and inserted codes marked by _CYL
2. DeltaPlasticStrain.c, same as above
3. MobilityLaw_BCC_0.c, same as above
4. MobilityLaw_FCC_0.c, same as above
5. InputSanity.c, same as above
6. LocalSegForces.c, same as above
7. NodeForce.c, same as above
8. PartialForces.c, same as above
9. Plot.c, same as above
10. Remesh.c, same as above
11. Mobility.h, removed from Include/
12. AllSegmentStress.c, copied from thinfilm/ and inserted codes marked by _CYL
13. CrossSlipFCC.c, CrossSlipBCC.c, WriteAtomEye.c, added to makefile
14. Include/CYL.h, added declaration of function Init3x3()
15. Include/ParadisProto.h, copied from thinfilm/

To do list:
1. remove global pointer cylinder0, modify Include/ParadisProto.h
2. create similar functions as those defined in Include/TF.h, TF_Util.c
3. do some test in the JMPS/image stress paper (ILL Ryu?)
4. merge final_matrices.c, greenstress_a.c, gridstress.c, m2invmatrix.c, etc
to CYL_Util.c
5. remove SegmentStress.c, replace by StressDueToSeg() in NodeForce.c
6. modify DeltaPlasticStrain.c ?
7. other fixes mentioned in README.txt (Sylvie)
8. check if stress field on surface is calculated corrected in PBC, similar to thinfilm/
9. double check NodeForce(), Yoffe stress, virtual segments, compare with thinfilm/
10. double check surface mobility law
11. double check surface remesh rule

----------------------------------------------------------------------------------
Sylvie.  Fri, May 8 2009

- I have not updated MobilityLaw_BCC_0a.c. 

- In NodeVelocity, only the new relax mobility law was added. It is now in the scr code, so
this file should not be changed anymore. Is this correct?

- Same for DeltaPlasticStrain.c

- Created a new variable called LenVirtualSeg. It was hard coded before.

- ComputeCYLSegSigbRem has been changed. Check it.

- Check Bending.

- I have not removed the cylinder0, home0 for now

- check virtual_segment_force: it was modified to get PARTIAL...
I have OrderNodes and not NodeOwnsSeg?????


- I have replaced hard array by dynamic ones. Check the order in which Tr, Tq and Tz are filled in in cylinder.c.

