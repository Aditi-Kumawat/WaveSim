clear; clc;
alpha_vect = [1.9 4 5.5 6.3 6.8 7.8];
beta_vect = [1.0 2.0 3.2 3.6 3.9 4.5];
rho_vect = [2.1 2.4 2.7 2.8 2.9 3.3];
h_vect = [0.5 1.0 2.5 23.0 13.0];
d_J=2.75;
[alpha_vect,beta_vect,rho_vect,h_vect,layer2split]=...
    fns_RmatGen.splitLayer(alpha_vect,beta_vect,rho_vect,h_vect,d_J);
disp(h_vect);
mu_vect = rho_vect .* beta_vect.^2;
%% k_vect range (divided in N_laps)
dk = 1.25e-4;
kmin=0;
kmax=15;
N_laps =960;
[k_mat_final, k_vect_last] = fns_RmatGen.setupKMat1(dk,kmin,kmax,N_laps);
%% Frequency range
f_hz_vect = 0.01:0.01:5;
omega_vect = 2 * pi * f_hz_vect;
o_I = 1e-3;
%%
F_x=1;
F_y=1;
F_z=1;
m_vect=[1 -1 0];

%% saving Rmat in this directory
dir_path = 'save data/R_mat_dJ2pt75';
if ~exist(dir_path, 'dir')
    mkdir(dir_path);
end
%% R matrix and r values computation
for i_m=1:length(m_vect)
    m=m_vect(i_m);
    [F_S,F_T,F_R]=fns_RmatGen.evalF_STR(m,F_x,F_y,F_z);
    for i_laps = 1:N_laps-1
        tic;
        fns_RmatGen.processLapwithsource(k_mat_final(i_laps, :),omega_vect,...
            o_I,i_laps,alpha_vect, beta_vect, h_vect, mu_vect,...
            dir_path,layer2split,F_S,F_R,F_T,m);
        display(strcat('lap=',num2str(i_laps)));

        elapsedTime = toc;
        fprintf('lap processing time: %f seconds.\n', elapsedTime);
    end
    %% R matrix and r values computation for last chunk
    fns_RmatGen.processLapwithsource(k_vect_last, omega_vect, o_I,N_laps,...
        alpha_vect, beta_vect, h_vect, mu_vect,dir_path,layer2split,F_S,F_R,F_T,m);
    display(strcat('lap=last'));
end





