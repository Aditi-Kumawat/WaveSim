clear;clc;close all
%%
r=10;
%% k_vect range (divided in N_laps)
dk = 1.25e-4;
kmin=1.25e-4;
kmax=20;
N_laps =2000;
k_laps= fns_RmatGen.setupKMat1(dk,kmin,kmax,N_laps);
kmax_new = k_laps(end,end);
%%
dir = 'save data/fSRT_Hmg_dJ_1e-5';
%% Frequency range
f_vect = 0.01:0.01:5;
omega_vect = 2 * pi * f_vect;
o_I = 1e-3;
tic;
%%
tic
for i_omg=1:length(omega_vect)
    omega=omega_vect(i_omg);
    f=f_vect(i_omg);
    display(strcat('f=',num2str(f),'Hz'));

    [fSRmat, fRRmat] = fns_GreenFnGen.load_fSRT_Hmg(f,dir);
    % [fSR_last, fRR_last, ~, ~] = fns_GreenFnGen.load_fSRT_last_Hmg(f,dir);
    Iz1=0;
    Ir1=0;
    for il=1:N_laps-1
        k=k_laps(il,:);
        J0 = besselj(0,k*r);
        J1 = besselj(1,k*r);
        fSR=fSRmat(il,:);
        fRR=fRRmat(il,:);
        Izint=k.*fRR.*J0;
        Iz(i_omg)=(trapz(k,Izint))+Iz1;
        Irint=k.*fSR.*J1;
        Ir(i_omg)=(trapz(k,Irint))+Ir1;
        Iz1=Iz(i_omg);
        Ir1=Ir(i_omg);
    end
    % k=k_vect_last;
    % J0 = besselj(0,k*r);
    % J1 = besselj(1,k*r);
    % fSR=fSR_last;
    % fRR=fRR_last;
    % Izint=k.*fRR.*J0;
    % Iz(i_omg)=(trapz(k,Izint))+Iz1;
    % Irint=k.*fSR.*J1;
    % Ir(i_omg)=(trapz(k,Irint))+Ir1;
end
elapsed_time = toc;
fprintf('run time dk integration: %.4f seconds\n', elapsed_time);

%%
t_vect=0:0.01:20;
T=1;
A=1;
s_vect=pi*A*T*(1+exp(1i*omega_vect*T))./(pi^2-omega_vect.^2*T^2);
idx_s_vect=find(abs(omega_vect-pi/T)<1e-12);
s_vect(idx_s_vect)=1i*A*T/2;
figure
plot(omega_vect,s_vect)
for i=1:length(t_vect)
    t=t_vect(1,i);
    integrand_hmg=s_vect.*Iz.*exp(-1i*omega_vect*t)/2/pi;
    I_hmg=exp(o_I*t)/pi*real(trapz(omega_vect,integrand_hmg));
    u_zz_t_vect_hmg(:,i)=I_hmg.';
    display(strcat('t=',num2str(t),'s'));
end
figure;
plot(t_vect,u_zz_t_vect_hmg)
