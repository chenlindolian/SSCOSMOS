# SSCOSMOS
**SMV kernel and Laplacian kernel based SSCOSMOS**

(paper link)

Single-step calculation of susceptibility with multiple orientation sampling (SS-COSMOS) was developed to jointly perform background field removal and COSMOS dipole inversion, which could serve as a better standard than traditional multi-step COSMOS (MS-COSMOS) for gauging single-orientation single-step quantitative susceptibilty mapping (SSQSM) methods. 

Two types of SSCOSMOS are tested and availabe here, i.e. SMV-kernel-based SSCOSMOS (with varialbe spherical-mean-value kernels), and Laplacian-kernel-based SSCOSMOS. The input should be MRI phase data acquired with GRE sequence at three or more head orientations.

- **SScosmos_SMV_lsmr.m** conducts SMV kernel based SSCOSMOS, the input needs to be unwrapped phase
- **SScosmos_Lap_lsmr.m** conducts Laplacian kernel based SSCOSMOS, the input can be raw wrapped phase.
- **testData.m** gives example usage of the Laplacian and SMV kernel based SSCOSMOS with some simulated phase data.
