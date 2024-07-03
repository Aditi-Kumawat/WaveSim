%% created by Aditi Kumawat:07.11.23
classdef fns_inptMatPara
    methods (Static)
        %%
        function soil_medium = select_soil_medium()
            options = {'homogeneous_HS', 'Wald1996_LHS', 'Munich_LHS','Poing_LHS'};
            [choice, isOk] = listdlg('PromptString', 'Select a soil medium:', ...
                'SelectionMode', 'single', ...
                'ListString', options,'ListSize',[150,150]);

            if isOk
                soil_medium = options{choice};
            else
                soil_medium = ''; % or handle the case where no selection is made
            end
        end
        %%
        function [alpha_vect,beta_vect,h_vect,JL,mu_vect,N,...
                alpha_L1,beta_L1,mu_L1,h_L1]=get_parameters(soil_medium,d_J)
            if strcmp(soil_medium, 'homogeneous_HS')
                alpha_vect = [2.4 2.4];
                beta_vect = [1.0 1.0];
                rho_vect = [2.3 2.3];
                h_vect = 0.5;
            elseif strcmp(soil_medium, 'Wald1996_LHS')
                alpha_vect = [1.9 4 5.5 6.3 6.8 7.8];
                beta_vect = [1.0 2.0 3.2 3.6 3.9 4.5];
                rho_vect = [2.1 2.4 2.7 2.8 2.9 3.3];
                h_vect = [0.5 1.0 2.5 23.0 13.0];
            elseif strcmp(soil_medium, 'Munich_LHS')
                alpha_vect = [850 2500 3000 3600 4400 6100]*1e-3;
                beta_vect = [350 950 1500 1850 2250 3500]*1e-3;
                rho_vect = [2300 2350 2500 2650 2900 2900]*1e-3;
                h_vect = [30,620,900,1500,600,1900]*1e-3;
            elseif strcmp(soil_medium, 'Poing_LHS')
                alpha_vect = [1050 2300 2450 3300 3900 4825 6100]*1e-3;
                beta_vect = [350 700 920 1700 2000 2475 3526]*1e-3;
                rho_vect = [2300 2300 2350 2480 2650 2820 2900]*1e-3;
                h_vect = [30 30 440 1000 900 855 2260]*1e-3;
            end
            [alpha_vect,beta_vect,rho_vect,h_vect,JL]=...
                fns_RmatGen.splitLayer(alpha_vect,beta_vect,rho_vect,...
                h_vect,d_J);
            mu_vect = rho_vect .* beta_vect.^2;
            N = length(alpha_vect) - 2;
            % Layer 1 properties
            alpha_L1=alpha_vect(1);
            beta_L1=beta_vect(1);
            mu_L1=mu_vect(1);
            h_L1=h_vect(1);
           
        end
        %% Directory name according to the soil medium and source depth
        function dir=form_dir(soil_medium,d_J)
            baseDir = 'SAVE_DATA';
            fldr_soilmedium = soil_medium;
            fldr_dpth=['fSRT_dJ', strrep(num2str(d_J), '.', 'pt'), '_',soil_medium];
            cd(baseDir)
            if ~exist(fldr_soilmedium, 'dir')
                mkdir(fldr_soilmedium);
            end
            cd(fldr_soilmedium)
            if ~exist(fldr_dpth, 'dir')
                mkdir(fldr_dpth);
            end
            cd(fullfile('..', '..'))
            %%-------------------------------------------------------%%
            dir = fullfile(baseDir, fldr_soilmedium, fldr_dpth);
            %%-------------------------------------------------------%%
        end

        function dir=form_dir_Iparas(soil_medium,d_J)
            baseDir = 'SAVE_DATA';
            fldr_soilmedium = soil_medium;
            fldr_dpth=['I_paras', strrep(num2str(d_J), '.', 'pt'), '_',soil_medium];
            cd(baseDir)
            if ~exist(fldr_soilmedium, 'dir')
                mkdir(fldr_soilmedium);
            end
            cd(fldr_soilmedium)
            if ~exist(fldr_dpth, 'dir')
                mkdir(fldr_dpth);
            end
            cd(fullfile('..', '..'))
            %%-------------------------------------------------------%%
            dir = fullfile(baseDir, fldr_soilmedium, fldr_dpth);
            %%-------------------------------------------------------%%
        end
    end
end