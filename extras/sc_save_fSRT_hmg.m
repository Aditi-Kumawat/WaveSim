clear; clc;
[alpha_vect,beta_vect,h_vect,JL,mu_vect,N,...
    alpha_L1,beta_L1,mu_L1,h_L1]=fns_inptMatPara.para_hmg;
ha_inv_2by2 = @fns_RmatGen.get_2by2_inv_3D_mat;
%% k_vect range (divided in N_laps)
dk = 1.25e-4;
kmin=1.25e-4;
kmax=20;
N_laps =2000;
k_laps= fns_RmatGen.setupKMat1(dk,kmin,kmax,N_laps);
kmax_new = k_laps(end,end);
%% Frequency range
f_hz_vect = 0.01:0.01:5;
omega_vect = 2 * pi * f_hz_vect;
o_I = 1e-3;
%%
F_x=0;
F_y=0;
F_z=1;
m_vect=[0];

%% saving Rmat in this directory
dir = 'save data/fSRT_Hmg_dJ_1e-5';
if ~exist(dir, 'dir')
    mkdir(dir);
end
%% R matrix and r values computation
for i_f=1:length(omega_vect)
    omega=omega_vect(i_f);
    f=f_hz_vect(i_f);
    display(strcat('f=',num2str(f),'Hz'));
    [fSRmat, fRRmat, fSRprm_mat, fRRprm_mat, fSSmat, fRSmat,...
        fSSprm_mat, fRSprm_mat, fTTmat, fTTprm_mat] =...
        fns_save_fSRT1.initialize_matrices(k_laps);
    tic;
    for i_m=1:length(m_vect)
        m=m_vect(i_m);
        [F_S,F_T,F_R]=fns_RmatGen.evalF_STR(m,F_x,F_y,F_z);

        for il = 1:N_laps-1
            k_vect=k_laps(il,:);
            [Rmat,AuBu,AuBubar,r,Cu,Cubar]=...
                fns_save_fSRT1.getRecurPara(k_vect,omega,o_I,...
                alpha_vect, beta_vect, h_vect, mu_vect,JL,...
                F_S,F_R,F_T,N,ha_inv_2by2);
            if m == 0
                [fSRmat(il,:), fRRmat(il,:),fSRprm_mat(il,:),...
                    fRRprm_mat(il,:)]=fns_save_fSRT1.get_fSRTmat_m0(omega,o_I,k_vect,...
                    alpha_L1,beta_L1, mu_L1,h_L1,F_R,...
                    Rmat,AuBu,AuBubar,r,Cu,Cubar);
            elseif m == 1
                % elseif m == 1 || m == -1
                %check if the values are same for m=+-1
                [fSSmat(il,:), fRSmat(il,:),fSSprm_mat(il,:),...
                    fRSprm_mat(il,:),fTTmat(il,:),fTTprm_mat(il,:)]=...
                    fns_save_fSRT1.get_fSRTmat_m1(omega,o_I,k_vect,alpha_L1,...
                    beta_L1, mu_L1,h_L1,F_S,F_T,Rmat,...
                    AuBu,AuBubar,r,Cu,Cubar);
            end
            % figure
            % plot(k_vect,fRRmat(il,:))
        end
        % %% R matrix and r values computation for last chunk
        % [Rmat,AuBu,AuBubar,r,Cu,Cubar]=...
        %     fns_save_fSRT1.getRecurPara(k_vect_last,omega,o_I,...
        %     alpha_vect, beta_vect, h_vect, mu_vect,JL,...
        %     F_S,F_R,F_T,N,ha_inv_2by2);
        % if m == 0
        %     [fSR_last, fRR_last,fSRprm_last,fRRprm_last]=...
        %         fns_save_fSRT1.get_fSRTmat_m0(omega,o_I,k_vect_last,...
        %         alpha_L1,beta_L1, mu_L1,h_L1,F_R,...
        %         Rmat,AuBu,AuBubar,r,Cu,Cubar);
        % elseif m == 1
        %     % elseif m == 1 || m == -1
        %     %check if the values are same for m=+-1
        %     [fSS_last, fRS_last,fSSprm_last,fRSprm_last,...
        %         fTT_last,fTTprm_last]=fns_save_fSRT1.get_fSRTmat_m1(...
        %         omega,o_I,k_vect_last,alpha_L1,beta_L1, mu_L1,h_L1,F_S,F_T,...
        %         Rmat,AuBu,AuBubar,r,Cu,Cubar);
        % end
    end
    fns_save_fSRT1.save_fSRT_HS(fSRmat, fRRmat,f,dir)

    % fns_save_fSRT1.save_fSRT_last_HS(fSR_last,fRR_last,fSRprm_last,...
    %     fRRprm_last,f,dir)
    elapsedTime = toc;
    fprintf('processing time: %f seconds.\n', elapsedTime);
end



