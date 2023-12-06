clear; clc;
[alpha_vect,beta_vect,h_vect,JL,mu_vect,N,...
    alpha_L1,beta_L1,mu_L1,h_L1]=fns_inptMatPara.layer_para_Poing;
ha_inv_2by2 = @fns_RmatGen.get_2by2_inv_3D_mat;

%%
r=5;
%% Frequency range
f = 3;
omega = 2 * pi * f;
o_I = 1e-3;
%%
F_x=0;
F_y=0;
F_z=1;
m=0;
[F_S,F_T,F_R] = fns_RmatGen.evalF_STR(m,F_x,F_y,F_z);
%% Test#1
% k_vect range (divided in N_laps)
dk = 1.25e-4;
kmin = 1.25e-4;
kmax = 10;
N_laps = 1000;
[k_laps] = fns_RmatGen.setupKMat1(dk,kmin,kmax,N_laps);

[~, fRRmat, ~, ~, ~, ~, ~, ~, ~, ~] = fns_save_fSRT1.initialize_matrices(k_laps);

Iz_accum = zeros(1, N_laps);  
tic;
parfor il = 1:N_laps
    k_vect = k_laps(il,:);
    J0 = besselj(0, k_vect*r);
    [Rmat,AuBu,AuBubar,rmat,Cu,Cubar] = fns_save_fSRT1.getRecurPara(k_vect,omega,o_I,...
        alpha_vect, beta_vect, h_vect, mu_vect,JL, F_S,F_R,F_T,N,ha_inv_2by2);
    [~, fRRmat(il,:), ~, ~] = fns_save_fSRT1.get_fSRTmat_m0(omega,o_I,k_vect,...
        alpha_L1,beta_L1, mu_L1,h_L1,F_R, Rmat,AuBu,AuBubar,rmat,Cu,Cubar);
    fRR = fRRmat(il,:);
    Iz_accum(il) = trapz(k_vect, k_vect.*fRR.*J0);  
end

Iz_test1 = sum(Iz_accum);  

elapsedTime = toc;
fprintf('processing time: %f seconds.\n', elapsedTime);




%% Test#2
% k_vect range (divided in N_laps)
dk = 1.25e-5;
kmin = 1.25e-5;
kmax = 10;
N_laps = 1000;
[k_laps] = fns_RmatGen.setupKMat1(dk,kmin,kmax,N_laps);

[~, fRRmat, ~, ~, ~, ~, ~, ~, ~, ~] = fns_save_fSRT1.initialize_matrices(k_laps);

Iz_accum = zeros(1, N_laps);  
tic;
parfor il = 1:N_laps
    k_vect = k_laps(il,:);
    J0 = besselj(0, k_vect*r);
    [Rmat,AuBu,AuBubar,rmat,Cu,Cubar] = fns_save_fSRT1.getRecurPara(k_vect,omega,o_I,...
        alpha_vect, beta_vect, h_vect, mu_vect,JL, F_S,F_R,F_T,N,ha_inv_2by2);
    [~, fRRmat(il,:), ~, ~] = fns_save_fSRT1.get_fSRTmat_m0(omega,o_I,k_vect,...
        alpha_L1,beta_L1, mu_L1,h_L1,F_R, Rmat,AuBu,AuBubar,rmat,Cu,Cubar);
    fRR = fRRmat(il,:);
    Iz_accum(il) = trapz(k_vect, k_vect.*fRR.*J0); 
end

Iz_test2 = sum(Iz_accum); 

elapsedTime = toc;
fprintf('processing time: %f seconds.\n', elapsedTime);


%%

dIuIR=abs(Iz_test2-Iz_test1)*100./abs(Iz_test2);
dIuIR_max=max(dIuIR);
dIuIR_min=min(dIuIR);