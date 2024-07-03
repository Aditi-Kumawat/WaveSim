%%
classdef fns_save_fSRT1
    methods (Static)
        %%
        function [Rmat,AuBu,AuBubar,r,Cu,Cubar]=getRecurPara(k_vect_lap, omega, o_I,...
                alpha_vect, beta_vect, h_vect, mu_vect,J,F_S,F_R,F_T,N,ha_inv_2by2)

            [Rmat,AuBu,AuBubar] = fns_get_recurPara.computeRplusAuBuMat(...
                omega, k_vect_lap, o_I,alpha_vect, beta_vect, h_vect,...
                mu_vect,N, ha_inv_2by2,J,F_S,F_R);
            [r,Cu,Cubar]=fns_get_recurPara.compute_rCu(omega,...
                k_vect_lap, o_I,beta_vect, h_vect, mu_vect,N,J,F_T);

        end
        %%
        function [fSRmat, fRRmat, fSRprm_mat, fRRprm_mat, ...
                fSSmat, fRSmat, fSSprm_mat, fRSprm_mat, ...
                fTTmat, fTTprm_mat] = initialize_matrices(k_mat_final)

            % Get the size of k_mat_final
            [M, N] = size(k_mat_final);

            % Initialize the matrices to zeros
            fSRmat = zeros(M, N);
            fRRmat = zeros(M, N);
            fSRprm_mat = zeros(M, N);
            fRRprm_mat = zeros(M, N);

            fSSmat = zeros(M, N);
            fRSmat = zeros(M, N);
            fSSprm_mat = zeros(M, N);
            fRSprm_mat = zeros(M, N);
            fTTmat = zeros(M, N);
            fTTprm_mat = zeros(M, N);
        end

        %%
        function [fSRmat, fRRmat,fSRprm_mat,fRRprm_mat]=get_fSRTmat_m0(...
                omega,o_I,k_vect,alpha_L1,beta_L1, mu_L1,...
                h_L1,F_R,RL13mat, AuBuL13vect,AuBubarL13vect,...
                rL1mat, CuL1vect,CuL1barvect)

            [u1S,u1R,u1Sprm,u1Rprm,~,~]=fns_save_fSRT1.getUL1(RL13mat,...
                AuBuL13vect,AuBubarL13vect, rL1mat, CuL1vect,CuL1barvect,alpha_L1,...
                beta_L1, mu_L1,omega, o_I, k_vect,h_L1);

            [fSRmat,fRRmat,fSRprm_mat,fRRprm_mat]=...
                fns_save_fSRT1.get_fSRT_m0(k_vect,F_R,u1S,u1R,u1Sprm,u1Rprm);
            % Check and report NaN values
            if any(isnan(fSRmat))
                disp('NaN detected in fSRmat.');
            end
            if any(isnan(fRRmat))
                disp('NaN detected in fRRmat.');
            end
            if any(isnan(fSRprm_mat))
                disp('NaN detected in fSRprm_mat.');
            end
            if any(isnan(fRRprm_mat))
                disp('NaN detected in fRRprm_mat.');
            end

        end
        %%
        function [fSSmat, fRSmat,fSSprm_mat,fRSprm_mat,fTTmat,...
                fTTprm_mat]=get_fSRTmat_m1(omega,o_I,k_vect,alpha_L1,beta_L1,...
                mu_L1,h_L1,F_S,F_T,RL13mat,...
                AuBuL13vect,AuBubarL13vect,...
                rL1mat, CuL1vect,CuL1barvect)

            [u1S,u1R,u1Sprm,u1Rprm,u1T,u1Tprm]=fns_save_fSRT1.getUL1(RL13mat,...
                AuBuL13vect,AuBubarL13vect, rL1mat, CuL1vect,CuL1barvect,alpha_L1,...
                beta_L1, mu_L1,omega, o_I, k_vect,h_L1);
            [fSSmat,fRSmat,fSSprm_mat,fRSprm_mat,...
                fTTmat,fTTprm_mat]=fns_save_fSRT1.get_fSRT_m1(k_vect,F_S,...
                F_T,u1S,u1R,u1Sprm,u1Rprm,u1T,u1Tprm);
            % Check and report NaN values
            if any(isnan(fSSmat))
                disp('NaN detected in fSRmat.');
            end
            if any(isnan(fRSmat))
                disp('NaN detected in fRRmat.');
            end
            if any(isnan(fSSprm_mat))
                disp('NaN detected in fSRprm_mat.');
            end
            if any(isnan(fRSprm_mat))
                disp('NaN detected in fRRprm_mat.');
            end
        end
        %%
        function [u1S,u1R,u1Sprm,u1Rprm,u1T,u1Tprm]=getUL1(RL13mat,...
                AuBuL13vect,AuBubarL13vect, rL1mat, CuL1vect,...
                CuL1barvect,alpha_L1,beta_L1, mu_L1,...
                omega, o_I, k_vect,h_L1)
            r11 = reshape(RL13mat(1,1,:), 1, []);
            r12 = reshape(RL13mat(1,2,:), 1, []);
            r21 = reshape(RL13mat(2,1,:), 1, []);
            r22 = reshape(RL13mat(2,2,:), 1, []);
            Au  = reshape(AuBuL13vect(1,1,:), 1, []);
            Bu  = reshape(AuBuL13vect(2,1,:), 1, []);
            Aubar = reshape(AuBubarL13vect(1,1,:), 1, []);
            Bubar = reshape(AuBubarL13vect(2,1,:), 1, []);
            %-------------------------------------------%
            r01=rL1mat;
            Cu=CuL1vect;
            Cubar=CuL1barvect;
            %-------------------------------------------%
            [p11,p13,p22,p24,p12,p14,p21,p23,p31,p33,p32,p34,p41,...
                p43,p42,p44,e11,e22]=fns_save_fSRT1.getlayer1props(...
                alpha_L1,beta_L1, mu_L1, omega, o_I, k_vect,...
                h_L1);
            [c11,c12,c21,c22]=fns_save_fSRT1.getcoeffABu1(p11,p13,p22,...
                p24,p12,p14,p21,p23,p31,p33,p32,p34,p41,p43,p42,p44,...
                e11,e22,r11,r12,r21,r22);
            u1S =c11.*Au+c12.*Bu;
            u1R =c21.*Au+c22.*Bu;
            u1Sprm =c11.*Aubar+c12.*Bubar;
            u1Rprm =c21.*Aubar+c22.*Bubar;
            u1T =2*(e22.*k_vect)./(1-e22.*r01).*Cu;
            u1Tprm =2*(e22.*k_vect)./(1-e22.*r01).*Cubar;
        end

        %%
        function [p11,p13,p22,p24,p12,p14,p21,p23,p31,p33,p32,p34,...
                p41,p43,p42,p44,e11,e22]=...
                getlayer1props(alpha, beta, mu, omega, o_I, k_vect, h)
            kBetaVect = (omega + 1i*o_I) ./ beta;
            kAlphaVect = (omega + 1i*o_I) ./ alpha;
            gammaVect = (kBetaVect.^2 - k_vect.^2).^(1/2);
            nuVect = (kAlphaVect.^2 - k_vect.^2).^(1/2);
            chiVect = k_vect.^2 - gammaVect.^2;

            p11 = k_vect;
            p13 = k_vect;
            p22 = k_vect;
            p24 = k_vect;

            p12 = 1i.*gammaVect;
            p14 = -1i.*gammaVect;
            p21 = 1i.*nuVect;
            p23 = -1i.*nuVect;

            p31 = 2*1i.*mu.*k_vect.*nuVect;
            p33 = -2*1i.*mu.*k_vect.*nuVect;
            p32 = mu*chiVect;
            p34 = mu*chiVect;

            p41 = mu*chiVect;
            p43 = mu*chiVect;
            p42 = 2*1i.*mu.*k_vect.*gammaVect;
            p44 = -2*1i.*mu.*k_vect.*gammaVect;

            %%
            e11=exp(1i.*nuVect*h);
            e22=exp(1i.*gammaVect*h);
        end
        %%
        function [c11,c12,c21,c22]=getcoeffABu1(p11,p13,p22,p24,p12,p14,p21,p23,p31,p33,p32,p34,p41,...
                p43,p42,p44,e11,e22,r11,r12,r21,r22)
            c11=e11.*(((p14 + p44.*((p11.*p32)./(p31.*p42 - p32.*p41) - (p12.*p31)./(p31.*p42 - p32.*p41)) - p34.*((p11.*p42)./(p31.*p42 - p32.*p41) - (p12.*p41)./(p31.*p42 - p32.*p41))).*(p31.*p43 - p33.*p41 + e22.*p31.*p42.*r21 - e22.*p32.*p41.*r21))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21) - ((p13 + p43.*((p11.*p32)./(p31.*p42 - p32.*p41) - (p12.*p31)./(p31.*p42 - p32.*p41)) - p33.*((p11.*p42)./(p31.*p42 - p32.*p41) - (p12.*p41)./(p31.*p42 - p32.*p41))).*(p31.*p44 - p34.*p41 + e22.*p31.*p42.*r22 - e22.*p32.*p41.*r22))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21));
            c12=e22.*(((p14 + p44.*((p11.*p32)./(p31.*p42 - p32.*p41) - (p12.*p31)./(p31.*p42 - p32.*p41)) - p34.*((p11.*p42)./(p31.*p42 - p32.*p41) - (p12.*p41)./(p31.*p42 - p32.*p41))).*(p32.*p43 - p33.*p42 - e11.*p31.*p42.*r11 + e11.*p32.*p41.*r11))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21) - ((p13 + p43.*((p11.*p32)./(p31.*p42 - p32.*p41) - (p12.*p31)./(p31.*p42 - p32.*p41)) - p33.*((p11.*p42)./(p31.*p42 - p32.*p41) - (p12.*p41)./(p31.*p42 - p32.*p41))).*(p32.*p44 - p34.*p42 - e11.*p31.*p42.*r12 + e11.*p32.*p41.*r12))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21));
            c21=e11.*(((p24 + p44.*((p21.*p32)./(p31.*p42 - p32.*p41) - (p22.*p31)./(p31.*p42 - p32.*p41)) - p34.*((p21.*p42)./(p31.*p42 - p32.*p41) - (p22.*p41)./(p31.*p42 - p32.*p41))).*(p31.*p43 - p33.*p41 + e22.*p31.*p42.*r21 - e22.*p32.*p41.*r21))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21) - ((p23 + p43.*((p21.*p32)./(p31.*p42 - p32.*p41) - (p22.*p31)./(p31.*p42 - p32.*p41)) - p33.*((p21.*p42)./(p31.*p42 - p32.*p41) - (p22.*p41)./(p31.*p42 - p32.*p41))).*(p31.*p44 - p34.*p41 + e22.*p31.*p42.*r22 - e22.*p32.*p41.*r22))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21));
            c22=e22.*(((p24 + p44.*((p21.*p32)./(p31.*p42 - p32.*p41) - (p22.*p31)./(p31.*p42 - p32.*p41)) - p34.*((p21.*p42)./(p31.*p42 - p32.*p41) - (p22.*p41)./(p31.*p42 - p32.*p41))).*(p32.*p43 - p33.*p42 - e11.*p31.*p42.*r11 + e11.*p32.*p41.*r11))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21) - ((p23 + p43.*((p21.*p32)./(p31.*p42 - p32.*p41) - (p22.*p31)./(p31.*p42 - p32.*p41)) - p33.*((p21.*p42)./(p31.*p42 - p32.*p41) - (p22.*p41)./(p31.*p42 - p32.*p41))).*(p32.*p44 - p34.*p42 - e11.*p31.*p42.*r12 + e11.*p32.*p41.*r12))./(p33.*p44 - p34.*p43 - e11.*p31.*p43.*r12 + e11.*p31.*p44.*r11 + e11.*p33.*p41.*r12 - e11.*p34.*p41.*r11 - e22.*p32.*p43.*r22 + e22.*p32.*p44.*r21 + e22.*p33.*p42.*r22 - e22.*p34.*p42.*r21 + e11.*e22.*p31.*p42.*r11.*r22 - e11.*e22.*p31.*p42.*r12.*r21 - e11.*e22.*p32.*p41.*r11.*r22 + e11.*e22.*p32.*p41.*r12.*r21));
        end
        %%
        function [fSR,fRR,fSRprm,fRRprm] = ...
                get_fSRT_m0(k_vect,F_R,u1S,u1R,u1Sprm,u1Rprm)
            fSR = u1S ./ k_vect ./ F_R;
            fRR = u1R ./ k_vect ./ F_R;
            fSRprm = u1Sprm ./ k_vect ./ F_R;
            fRRprm = u1Rprm ./ k_vect ./ F_R;
        end

        %%
        function [fSS,fRS,fSSprm,fRSprm,fTT,fTTprm] = ...
                get_fSRT_m1(k_vect,F_S,F_T,u1S,u1R,u1Sprm,u1Rprm,u1T,u1Tprm)
            fSS = u1S ./ k_vect ./ F_S;
            fRS = u1R ./ k_vect ./ F_S;
            fSSprm = u1Sprm ./ k_vect ./ F_S;
            fRSprm = u1Rprm ./ k_vect ./ F_S;
            fTT = u1T ./ k_vect ./ F_T;
            fTTprm = u1Tprm ./ k_vect ./ F_T;
        end
        %%
        function save_fSRT(fSRmat, fRRmat,fSRprm_mat,fRRprm_mat,fSSmat,...
                fRSmat,fSSprm_mat,fRSprm_mat,fTTmat,fTTprm_mat,f,dir)
            %% saving f_SRT

            filename1 = fullfile(dir, sprintf('fSR_f_%d_Hz.mat',f));
            save(filename1, 'fSRmat');

            filename2 = fullfile(dir, sprintf('fRR_f_%d_Hz.mat',f));
            save(filename2, 'fRRmat');

            filename3 = fullfile(dir, sprintf('fSRprm_f_%d_Hz.mat',f));
            save(filename3, 'fSRprm_mat');

            filename4 = fullfile(dir, sprintf('fRRprm_f_%d_Hz.mat',f));
            save(filename4, 'fRRprm_mat');

            filename5 = fullfile(dir, sprintf('fSS_f_%d_Hz.mat',f));
            save(filename5, 'fSSmat');

            filename6 = fullfile(dir, sprintf('fRS_f_%d_Hz.mat',f));
            save(filename6, 'fRSmat');

            filename7 = fullfile(dir, sprintf('fSSprm_f_%d_Hz.mat',f));
            save(filename7, 'fSSprm_mat');

            filename8 = fullfile(dir, sprintf('fRSprm_f_%d_Hz.mat',f));
            save(filename8, 'fRSprm_mat');

            filename9 = fullfile(dir, sprintf('fTT_f_%d_Hz.mat',f));
            save(filename9, 'fTTmat');

            filename10 = fullfile(dir, sprintf('fTTprm_f_%d_Hz.mat',f));
            save(filename10, 'fTTprm_mat');
        end
        %%
        function save_fSRT_HS(fSRmat, fRRmat,f,dir)
            %% saving f_SRT

            filename1 = fullfile(dir, sprintf('fSR_f_%d_Hz.mat',f));
            save(filename1, 'fSRmat');

            filename2 = fullfile(dir, sprintf('fRR_f_%d_Hz.mat',f));
            save(filename2, 'fRRmat')

        end
       
    end
end