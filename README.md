# Source_Wavefield_Simulation
Generation of synthetic ground motions for a 1D soil profile through a fast and simplified approach. The approach idealized the soil as an N+1 layered soil medium with an underlying half-space. Each layer is assumed to be elastic, homogeneous, and isotropic, with constant material properties across the layer. The layers are assumed to be free of the body forces. The source has been modelled as a strike-slip fault and is located in the jth layer. The P, SV, and SH amplitudes are defined for each layer, and then the boundary conditions (radiation at infinite depth, displacement continuity across layers, traction continuity across layers, free-surface) are used to relate those amplitudes and derive the recurrence relation and the Green's function. The analysis used is described in more detail in the papers by Singla and Gupta 2019 (https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/109/1/433/568060/Surface-Rotations-Due-to-Kinematic-Shear?redirectedFrom=fulltext) and Pei et al. 2008 (https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/98/1/280/341880/Improvements-on-Computation-of-Phase-Velocities-of?redirectedFrom=fulltext). 
![image](https://github.com/Aditi-Kumawat/Source_Wavefield_Simulation/assets/72736535/4b977241-f0f9-44fe-8aeb-b3d0c3883af6)

Code Running Sequence:
1. sc_save_fSRT_ll.m: generates and saves soil medium+source-depth dependent parameters fSS; fSR; fRS; fRR, and fTT (eqn 41) and f′SS; f′SR; f′RS; f′RR, and f′TT(appendix), check for soil medium, depth, k_vect (N_laps)
2. sc_save_I_paras.m: saves I_r^(1)(r,\omega), I_r^(2)(r,\omega), I_r^(3)(r,\omega), I_phi^(1)(r,\omega), I_phi^(2)(r,\omega), I_z^(1)(r,\omega), I_z^(2)(r,\omega), I_z^(3)(r,\omega) (eqns: 56-63)
3. sc_velCYLcoordi.m: generates ground motion displacement using Eqns (46, 48-51) at given epicentral distances and source properties.

*Equation numbers correspond to the study
