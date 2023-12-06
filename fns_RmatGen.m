%%
classdef fns_RmatGen
    methods (Static)
        %%
        function [alpha_vect,beta_vect,rho_vect,h_vect,layer2split]=...
                splitLayer(alpha_vect,beta_vect,rho_vect,h_vect,d_J)
            % If d_J is 0, then simply return without splitting
            if d_J == 0
                layer2split='none';
                return;
            end
            cumulative_h = 0;  % initialize cumulative depth
            layer2split = 0;  % initialize layer index
            % Find the layer to split
            for i = 1:length(h_vect)
                cumulative_h = cumulative_h + h_vect(i);
                if cumulative_h >= d_J
                    layer2split = i;
                    break;
                end
            end
            % If d_J is within the existing layers
            if layer2split > 0
                % Calculate new thickness for split layers
                h_before_dJ = d_J - (cumulative_h - h_vect(layer2split));
                h_after_dJ = h_vect(layer2split) - h_before_dJ;
                % Adjust the h_vect, alpha_vect, beta_vect, and rho_vect
                h_vect = [h_vect(1:layer2split-1) h_before_dJ...
                    h_after_dJ h_vect(layer2split+1:end)];
                alpha_vect = [alpha_vect(1:layer2split)...
                    alpha_vect(layer2split) alpha_vect(layer2split+1:end)];
                beta_vect = [beta_vect(1:layer2split)...
                    beta_vect(layer2split) beta_vect(layer2split+1:end)];
                rho_vect = [rho_vect(1:layer2split)...
                    rho_vect(layer2split) rho_vect(layer2split+1:end)];
                % If d_J is larger than the total thickness
            else
                extra_thickness = d_J - cumulative_h;
                layer2split=length(alpha_vect);
                % Append the new layer to h_vect
                h_vect = [h_vect, extra_thickness];
                % Add the properties of the medium underneath
                alpha_vect = [alpha_vect, alpha_vect(end)];
                beta_vect = [beta_vect, beta_vect(end)];
                rho_vect = [rho_vect, rho_vect(end)];
            end
        end
        %%
         %%
        function [F_S,F_T,F_R]=evalF_STR(m,F_x,F_y,F_z)
            switch m
                case 1
                    F_T=-(1/4/pi)*(1i*F_x+F_y);
                    F_S=(1/4/pi)*(F_x-1i*F_y);
                    F_R=0;
                case -1
                    F_T=-(1/4/pi)*(1i*F_x-F_y);
                    F_S=-(1/4/pi)*(F_x+1i*F_y);
                    F_R=0;
                case 0
                    F_T=0;
                    F_S=0;
                    F_R=F_z*(1/2/pi);
            end
        end
        %%
        function k_mat = setupKMat1(dk,kmin,kmax,N_laps)
            k_vect = kmin:dk:kmax;
            lap_length = (length(k_vect) / N_laps);
            if mod(lap_length, 1) ~= 0
                error('Lap length is not an integer. Exiting the function.');
            end
            k_mat = zeros(N_laps, lap_length);

            for i_lap = 1:N_laps
                start_idx = (i_lap - 1) * lap_length + 1;
                end_idx = i_lap * lap_length;
                k_mat(i_lap, :) = k_vect(start_idx:end_idx);
            end
        end

        %%
        function processLapwithsource(k_vect_lap, omega_vect, o_I, lap_num,...
                alpha_vect, beta_vect, h_vect, mu_vect,dir_path,J,F_S,F_R,F_T,m)

            N = length(alpha_vect) - 2;
            ha_inv_2by2 = @fns_RmatGen.get_2by2_inv_3D_mat;
            Rmat_layer1 = zeros(length(omega_vect),length(k_vect_lap), 4);
            AuBu_layer1 = zeros(length(omega_vect),length(k_vect_lap), 2);
            AuBubar_layer1 = zeros(length(omega_vect),length(k_vect_lap), 2);
            r_layer1 = zeros(length(omega_vect),length(k_vect_lap));
            Cu_layer1 = zeros(length(omega_vect),length(k_vect_lap));
            Cubar_layer1 = zeros(length(omega_vect),length(k_vect_lap));
            for i_omg = 1:length(omega_vect)
                omega = omega_vect(i_omg);
                [Rmat,AuBu,AuBubar] = fns_get_recurPara.computeRplusAuBuMat(...
                    omega, k_vect_lap, o_I,alpha_vect, beta_vect, h_vect,...
                    mu_vect,N, ha_inv_2by2,J,F_S,F_R);
                [r,Cu,Cubar]=fns_get_recurPara.compute_rCu(omega,...
                    k_vect_lap, o_I,beta_vect, h_vect, mu_vect,N,J,F_T);

                Rmat_layer1(i_omg, :, 1) = Rmat(1, 1, :);
                Rmat_layer1(i_omg, :, 2) = Rmat(1, 2, :);
                Rmat_layer1(i_omg, :, 3) = Rmat(2, 1, :);
                Rmat_layer1(i_omg, :, 4) = Rmat(2, 2, :);
                AuBu_layer1(i_omg, :, 1) = AuBu(1, 1, :);
                AuBu_layer1(i_omg, :, 2) = AuBu(2, 1, :);
                AuBubar_layer1(i_omg, :, 1) = AuBubar(1, 1, :);
                AuBubar_layer1(i_omg, :, 2) = AuBubar(2, 1, :);
                r_layer1(i_omg,:)=r;
                Cu_layer1(i_omg,:)=Cu;
                Cubar_layer1(i_omg,:)=Cubar;
            end

            fns_RmatGen.save_recurpara(Rmat_layer1,AuBu_layer1,AuBubar_layer1,r_layer1,Cu_layer1,Cubar_layer1,lap_num,dir_path,m);
        end
        %%
        function save_recurpara(Rmat_layer1,AuBu_layer1,AuBubar_layer1,...
                r_layer1,Cu_layer1,Cubar_layer1,lap,dir_path,m)
            filename2 = fullfile(dir_path, sprintf('r_lap%d_m%d.mat',lap,m));
            save(filename2, 'r_layer1');
            filename1 = fullfile(dir_path, sprintf('AuBuvect_lap%d_m%d.mat',lap,m));
            save(filename1, 'AuBu_layer1');
            filename4 = fullfile(dir_path, sprintf('AuBubarvect_lap%d_m%d.mat',lap,m));
            save(filename4, 'AuBubar_layer1');
            filename = fullfile(dir_path, sprintf('Rmat_lap%d_m%d.mat',lap,m));
            save(filename, 'Rmat_layer1');
            filename3 = fullfile(dir_path, sprintf('Cu_lap%d_m%d.mat',lap,m));
            save(filename3, 'Cu_layer1');
            filename3 = fullfile(dir_path, sprintf('Cubar_lap%d_m%d.mat',lap,m));
            save(filename3, 'Cubar_layer1');
        end

        %% Evaluates the inverse of each slice of a 3D matrix
        function inv_3D=get_2by2_inv_3D_mat(mat_3D,length_3D)

            mat_3D_11=reshape(mat_3D(1,1,:),[1,length_3D]);
            mat_3D_12=reshape(mat_3D(1,2,:),[1,length_3D]);
            mat_3D_21=reshape(mat_3D(2,1,:),[1,length_3D]);
            mat_3D_22=reshape(mat_3D(2,2,:),[1,length_3D]);
            det_mat_3D=mat_3D_11.*mat_3D_22-mat_3D_12.*mat_3D_21;

            inv_3D(1,1,:)=mat_3D_22./det_mat_3D;
            inv_3D(1,2,:)=-mat_3D_12./det_mat_3D;
            inv_3D(2,1,:)=-mat_3D_21./det_mat_3D;
            inv_3D(2,2,:)=mat_3D_11./det_mat_3D;
        end
        %% INPUT
        %   M  : 3D array (m x n x p)
        %   RHS: 3D array (m x q x p)
        % OUTPUT
        %   X  : 3D array (n x q x p)
        %
        % Solve p systems of linear equations:
        %   M(:,:,k) * X(:,:,k) = RHS(:,:,k) for all k=1,2,...,p
        %
        % NOTE 1 -- Special call with squeezed RHS: common RHS for all systems.
        %   Input argument RHS can be squeezed as (m x q),
        %       M   is 3D array (m x n x p)
        %       RHS is 2D array (m x q)
        %       X   is 2D array (n x q x p)
        %       M(:,:,k) * X(:,:,k) = RHS for all k=1,2,...,p
        function X = SliceMultiSolver(M, RHS)
            RHS = permute(RHS, [1 3 2]);
            X = MultiSolver(M, RHS);
            if size(M,3)>1
                X = permute(X, [1 3 2]);
            end
        end
        %% INPUT
        %   M  : 3D array (m x n x p)
        %   X  : 3D array (n x q x p)
        % OUTPUT
        %   RHS: 3D array (m x q x p)
        %
        % Compute p matrix products:
        %   M(:,:,k) * X(:,:,k) = RHS(:,:,k) for all k=1,2,...,p
        %
        % NOTE -- Special call with X: common right matrix
        %   Input argument X can be contracted as (n x q), and automatically
        %   expanded by the function
        %       M   is 3D array (m x n x p)
        %       X   is 2D array (n x q)
        %       RHS is 3D array (m x q x p)
        %       M(:,:,k) * X = RHS(:,:,k) for all k=1,2,...,p
        function RHS = SliceMultiProd(M, X)
            X = permute(X, [1 3 2]);
            RHS = MultiProd(M, X);
            if size(M,3)>1
                RHS = permute(RHS, [1 3 2]);
            end
        end
    end
end