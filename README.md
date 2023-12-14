# Source_Wavefield_Simulation
Generation of synthetic ground motions for a 1D soil profile:
Code Running Sequence:
1. sc_save_fSRT_ll.m: check for soil medium, depth, k_vect (N_laps): generates and saves soil medium+source-depth dependent parameters fSS; fSR; fRS; fRR, and fTT (eqn 41) and f′SS; f′SR; f′RS; f′RR, and f′TT(appendix)
2. sc_save_I_paras.m: saves I_r^(1)(r,\omega), I_r^(2)(r,\omega), I_r^(3)(r,\omega), I_phi^(1)(r,\omega), I_phi^(2)(r,\omega), I_z^(1)(r,\omega), I_z^(2)(r,\omega), I_z^(3)(r,\omega) (eqns: 56-63)
3. sc_velCYLcoordi.m: generates ground motion displacement using eqns (46, 48-51) at given epicentral distances and source properties
