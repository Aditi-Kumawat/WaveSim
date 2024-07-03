%% created by Aditi Kumawat:07.11.23
classdef fns_GreenFnGen
    methods (Static)
        %%
        function [fSRmat, fRRmat, fSRprm_mat, fRRprm_mat, fSSmat, fRSmat,...
                fSSprm_mat, fRSprm_mat, fTTmat,fTTprm_mat] = load_fSRT(f,dir)
            filename1 = fullfile(dir, sprintf('fSR_f_%d_Hz.mat',f));
            data1 = load(filename1);
            fSRmat = data1.fSRmat;

            filename2 = fullfile(dir, sprintf('fRR_f_%d_Hz.mat',f));
            data2 = load(filename2);
            fRRmat = data2.fRRmat;

            filename3 = fullfile(dir, sprintf('fSRprm_f_%d_Hz.mat',f));
            data3 = load(filename3);
            fSRprm_mat = data3.fSRprm_mat;

            filename4 = fullfile(dir, sprintf('fRRprm_f_%d_Hz.mat',f));
            data4 = load(filename4);
            fRRprm_mat = data4.fRRprm_mat;

            filename5 = fullfile(dir, sprintf('fSS_f_%d_Hz.mat',f));
            data5 = load(filename5);
            fSSmat = data5.fSSmat;

            filename6 = fullfile(dir, sprintf('fRS_f_%d_Hz.mat',f));
            data6 = load(filename6);
            fRSmat = data6.fRSmat;

            filename7 = fullfile(dir, sprintf('fSSprm_f_%d_Hz.mat',f));
            data7 = load(filename7);
            fSSprm_mat = data7.fSSprm_mat;

            filename8 = fullfile(dir, sprintf('fRSprm_f_%d_Hz.mat',f));
            data8 = load(filename8);
            fRSprm_mat = data8.fRSprm_mat;

            filename9 = fullfile(dir, sprintf('fTT_f_%d_Hz.mat',f));
            data9 = load(filename9);
            fTTmat = data9.fTTmat;

            filename10 = fullfile(dir, sprintf('fTTprm_f_%d_Hz.mat',f));
            data10 = load(filename10);
            fTTprm_mat = data10.fTTprm_mat;
        end
        %%
        function [fSRmat, fRRmat] = load_fSRT_Hmg(f,dir)
            filename1 = fullfile(dir, sprintf('fSR_f_%d_Hz.mat',f));
            data1 = load(filename1);
            fSRmat = data1.fSRmat;

            filename2 = fullfile(dir, sprintf('fRR_f_%d_Hz.mat',f));
            data2 = load(filename2);
            fRRmat = data2.fRRmat;
        end
        %% saving I parameters  
        function save_I_r_phi_z(Ir1,Ir2,Ir3,Iph1,Iph2,Iz1,Iz2,Iz3,r,dir)

            filename1 = fullfile(dir, sprintf('Ir1_r_%d_m.mat',r));
            save(filename1, 'Ir1');

            filename2 = fullfile(dir, sprintf('Ir2_r_%d_m.mat',r));
            save(filename2, 'Ir2')

            filename3 = fullfile(dir, sprintf('Ir3_r_%d_m.mat',r));
            save(filename3, 'Ir3')

            filename4 = fullfile(dir, sprintf('Iph1_r_%d_m.mat',r));
            save(filename4, 'Iph1');

            filename5 = fullfile(dir, sprintf('Iph2_r_%d_m.mat',r));
            save(filename5, 'Iph2')

            filename6 = fullfile(dir, sprintf('Iz1_r_%d_m.mat',r));
            save(filename6, 'Iz1');

            filename7 = fullfile(dir, sprintf('Iz2_r_%d_m.mat',r));
            save(filename7, 'Iz2')

            filename8 = fullfile(dir, sprintf('Iz3_r_%d_m.mat',r));
            save(filename8, 'Iz3')

        end
%%
        function [Ir1, Ir2, Ir3, Iph1, Iph2, Iz1, Iz2, Iz3] = load_I_r_phi_z(r, dir)

            filename1 = fullfile(dir, sprintf('Ir1_r_%d_m.mat', r));
            temp = load(filename1);
            Ir1 = temp.Ir1;

            filename2 = fullfile(dir, sprintf('Ir2_r_%d_m.mat', r));
            temp = load(filename2);
            Ir2 = temp.Ir2;

            filename3 = fullfile(dir, sprintf('Ir3_r_%d_m.mat', r));
            temp = load(filename3);
            Ir3 = temp.Ir3;

            filename4 = fullfile(dir, sprintf('Iph1_r_%d_m.mat', r));
            temp = load(filename4);
            Iph1 = temp.Iph1;

            filename5 = fullfile(dir, sprintf('Iph2_r_%d_m.mat', r));
            temp = load(filename5);
            Iph2 = temp.Iph2;

            filename6 = fullfile(dir, sprintf('Iz1_r_%d_m.mat', r));
            temp = load(filename6);
            Iz1 = temp.Iz1;

            filename7 = fullfile(dir, sprintf('Iz2_r_%d_m.mat', r));
            temp = load(filename7);
            Iz2 = temp.Iz2;

            filename8 = fullfile(dir, sprintf('Iz3_r_%d_m.mat', r));
            temp = load(filename8);
            Iz3 = temp.Iz3;
        end

    end
end