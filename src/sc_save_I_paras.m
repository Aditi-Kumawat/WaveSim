% created by Aditi Kumawat:07.11.23
clear;clc;close all
%%
r_vect=[1 10];
%%-------------------------------------------------------%%
soil_medium = fns_inptMatPara.select_soil_medium();
disp(['selectedRfFldr: ', soil_medium])
%%-------------------------------------------------------%%
d_J=0.25; % Depth value
%%-------------------------------------------------------%%
dir=fns_inptMatPara.form_dir(soil_medium,d_J);
dir_Iparas=fns_inptMatPara.form_dir_Iparas(soil_medium,d_J);
%%-------------------------------------------------------%%
%% k_vect range (divided in N_laps)
filename_k = fullfile(dir_Iparas, 'k_laps.mat');
load(filename_k);
filename_Nlaps = fullfile(dir_Iparas, 'N_laps.mat');
load(filename_Nlaps);
%% Frequency range
filename_f = fullfile(dir_Iparas, 'f_vect.mat');
load(filename_f);

omega_vect = 2 * pi * f_vect;
o_I = 1e-3;
tic;

%% Start a parallel pool
if isempty(gcp('nocreate'))
    parpool;
end
%%
for i_r=1:length(r_vect)
    r=r_vect(i_r);
     display(strcat('r=',num2str(r),'m'));
    parfor i_omg=1:length(omega_vect)
        omega=omega_vect(i_omg);
        f=f_vect(i_omg);
        % display(strcat('f=',num2str(f),'Hz'));
               [fSRmat, fRRmat, fSRprm_mat, fRRprm_mat, fSSmat, fRSmat, fSSprm_mat,...
            fRSprm_mat, fTTmat,fTTprm_mat] = fns_GreenFnGen.load_fSRT(f,dir);
        Ir11=0;
        Ir21=0;
        Ir31=0;
        Iph11=0;
        Iph21=0;
        Iz11=0;
        Iz21=0;
        Iz31=0;
        for il=1:N_laps-1
            k=k_laps(il,:);
            J0 = besselj(0,k*r);
            J1 = besselj(1,k*r);
            J2 = besselj(2,k*r);
            J3 = besselj(3,k*r);

            fSR=fSRmat(il,:);
            fRR=fRRmat(il,:);
            fSRbar=fSRprm_mat(il,:);
            fRRbar=fRRprm_mat(il,:);
            fSS=fSSmat(il,:);
            fRS=fRSmat(il,:);
            fSSbar=fSSprm_mat(il,:);
            fRSbar=fRSprm_mat(il,:);
            fTT=fTTmat(il,:);
            fTTbar=fTTprm_mat(il,:);
            Ir1int=k.^2.*(fSS.*(J1-J3)+fTT.*(J1+J3));
            Ir2int=k.*((k.*fSR-fSSbar).*(J0-J2)-fTTbar.*(J0+J2));
            Ir3int=k.*(k.*fSS-2*fSRbar).*J1;
            Iph1int=k.^2.*(fSS.*(J1+J3)+fTT.*(J1-J3));
            Iph2int=k.*((fSSbar-k.*fSR).*(J0+J2)+fTTbar.*(J0-J2));
            Iz1int=k.^2.*fRS.*J2;
            Iz2int=k.*(k.*fRR-fRSbar).*J1;
            Iz3int=k.*(k.*fRS-2.*fRRbar).*J0;
            Ir1(i_omg)=(trapz(k,Ir1int))+Ir11;
            Ir2(i_omg)=(trapz(k,Ir2int))+Ir21;
            Ir3(i_omg)=(trapz(k,Ir3int))+Ir31;
            Iph1(i_omg)=(trapz(k,Iph1int))+Iph11;
            Iph2(i_omg)=(trapz(k,Iph2int))+Iph21;
            Iz1(i_omg)=(trapz(k,Iz1int))+Iz11;
            Iz2(i_omg)=(trapz(k,Iz2int))+Iz21;
            Iz3(i_omg)=(trapz(k,Iz3int))+Iz31;
            Ir11=Ir1(i_omg);
            Ir21=Ir2(i_omg);
            Ir31=Ir3(i_omg);
            Iph11=Iph1(i_omg);
            Iph21=Iph2(i_omg);
            Iz11=Iz1(i_omg);
            Iz21=Iz2(i_omg);
            Iz31=Iz3(i_omg);
        end
    end
    fns_GreenFnGen.save_I_r_phi_z(Ir1,Ir2,Ir3,Iph1,Iph2,Iz1,Iz2,Iz3,r,dir_Iparas)
end
elapsed_time = toc;
fprintf('run time: %.4f seconds\n', elapsed_time);