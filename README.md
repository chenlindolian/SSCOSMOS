# SSCOSMOS
**SMV kernel and Laplacian kernel based SSCOSMOS**

SSCOSMOS was developed for single step COSMOS (SSCOSMOS) to jointly perform background field removal and dipole inversion with multiple orientation sampling, which could serve as a better standard for gauging SSQSM methods. It is incorporated  multiple spherical mean value (SMV) kernels of various radii with the dipole inversion in SSCOSMOS. 

We also incorporated Laplacian kernel with the dipole inversion in SSCOSMOS, the phase unwrapping procedure can also be integrated into single-step QSM.


**SScosmos_Lap_lsmr.m** conduct Laplacian kernel based SSCOSMOS, the input can be raw wrapped phase.
**SScosmos_SMV_lsmr.m** conduct SMV kernel based SSCOSMOS, the input should be the unwrapped phase

The input data should be acquired at three or more head orientations which can be reconstructed for better performance than traditional multi-step COSMOS.

**testData.m** gives the example for Laplacian and SMV kernel based SSCOSMOS using simulation data.
