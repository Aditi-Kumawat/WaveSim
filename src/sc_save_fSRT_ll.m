% created by Aditi Kumawat:07.11.23
clear; clc;
TimeStart = tic;
%%-------------------------------------------------------%%
soil_medium = fns_inptMatPara.select_soil_medium();
disp(['selectedRfFldr: ', soil_medium])
%%-------------------------------------------------------%%
d_J=0.25; % Depth value
%%-------------------------------------------------------%%
[alpha_vect,beta_vect,h_vect,JL,mu_vect,N,...
    alpha_L1,beta_L1,mu_L1,h_L1]=fns_inptMatPara.get_parameters(soil_medium,d_J);
ha_inv_2by2 = @fns_RmatGen.get_2by2_inv_3D_mat;
%%-------------------------------------------------------%%
dir=fns_inptMatPara.form_dir(soil_medium,d_J); % saving parameters
dir_Iparas=fns_inptMatPara.form_dir_Iparas(soil_medium,d_J);
%%-------------------------------------------------------%%
%% k_vect range (divided in N_laps)
dk = 1.25e-4;
kmin=1.25e-4;
kmax=65;
N_laps =6500;
k_laps= fns_RmatGen.setupKMat1(dk,kmin,kmax,N_laps);
filename_k = fullfile(dir_Iparas, 'k_laps.mat');
save(filename_k, 'k_laps');
filename_Nlaps = fullfile(dir_Iparas, 'N_laps.mat');
save(filename_Nlaps, 'N_laps');
%% Frequency range
f_vect = 0.01:0.01:5;
omega_vect = 2 * pi * f_vect;
o_I = 1e-3;
filename_f = fullfile(dir_Iparas, 'f_vect.mat');
save(filename_f, 'f_vect');
%%
F_x=1;
F_y=1;
F_z=1;
m_vect=[1 0];

%% Start a parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

%% fSRT computation
parfor i_f=1:length(omega_vect)
    omega=omega_vect(i_f);
    f=f_vect(i_f);
    % display(strcat('f=',num2str(f),'Hz'));
    [fSRmat,fRRmat,fSRprm_mat,fRRprm_mat,fSSmat,fRSmat,fSSprm_mat,...
        fRSprm_mat, fTTmat, fTTprm_mat] =...
        fns_save_fSRT1.initialize_matrices(k_laps);
    for i_m=1:length(m_vect)
        m=m_vect(i_m);
        [F_S,F_T,F_R]=fns_RmatGen.evalF_STR(m,F_x,F_y,F_z);

        for il = 1:N_laps
            k_vect=k_laps(il,:);
            [Rmat,AuBu,AuBubar,r,Cu,Cubar]=...
                fns_save_fSRT1.getRecurPara(k_vect,omega,o_I,...
                alpha_vect, beta_vect, h_vect, mu_vect,JL,...
                F_S,F_R,F_T,N,ha_inv_2by2);
            if m == 0
                [fSRmat(il,:), fRRmat(il,:),fSRprm_mat(il,:),...
                    fRRprm_mat(il,:)]=fns_save_fSRT1.get_fSRTmat_m0(...
                    omega,o_I,k_vect,alpha_L1,beta_L1,mu_L1,h_L1,F_R,...
                    Rmat,AuBu,AuBubar,r,Cu,Cubar);
            elseif m == 1
                [fSSmat(il,:), fRSmat(il,:),fSSprm_mat(il,:),...
                    fRSprm_mat(il,:),fTTmat(il,:),fTTprm_mat(il,:)]=...
                    fns_save_fSRT1.get_fSRTmat_m1(omega,o_I,k_vect,...
                    alpha_L1,beta_L1, mu_L1,h_L1,F_S,F_T,Rmat,...
                    AuBu,AuBubar,r,Cu,Cubar);
            end
        end
    end

    fns_save_fSRT1.save_fSRT(fSRmat,fRRmat,fSRprm_mat,fRRprm_mat,...
        fSSmat,fRSmat,fSSprm_mat,fRSprm_mat,fTTmat,fTTprm_mat,f,dir)
end
totalElapsedTime = toc(TimeStart);
fprintf('Total processing time: %f seconds.\n', totalElapsedTime);
