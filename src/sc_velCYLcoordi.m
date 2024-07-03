% clear;clc;close all
%%
r=1;
phi_s=deg2rad(30);          %strike
delta_s=deg2rad(60);        %dip
lamda_s=deg2rad(50);       %slip or rake
M_w=4;
M_0=10^(1.5*(M_w+10.7));
% M_0=1e22;
%%-------------------------------------------------------%%
soil_medium = fns_inptMatPara.select_soil_medium();
disp(['selectedRfFldr: ', soil_medium])
%%-------------------------------------------------------%%
d_J=0.25; % Depth value
%%-------------------------------------------------------%%
dir=fns_inptMatPara.form_dir_Iparas(soil_medium,d_J);
%%-------------------------------------------------------%%
%% Frequency range
filename_f = fullfile(dir, 'f_vect.mat');
load(filename_f);
omega_vect = 2 * pi * f_vect;
o_I = 1e-3;
%%
[Ir1, Ir2, Ir3, Iph1, Iph2, Iz1, Iz2, Iz3] = fns_GreenFnGen.load_I_r_phi_z(r, dir);

gu1_r=sin(2*phi_s)/(4*pi).*(Ir1);
gu1_ph=cos(2*phi_s)/(4*pi).*(Iph1);
gu1_z=sin(2*phi_s)/(2*pi).*(Iz1);

gu2_r=cos(phi_s)/(4*pi).*(Ir2);
gu2_ph=sin(phi_s)/(4*pi).*(Iph2);
gu2_z=cos(phi_s)/(2*pi).*(Iz2);

gu3_r=(1/(4*pi)).*Ir3 - (cos(2*phi_s)/(8*pi))*Ir1;
gu3_ph=(sin(2*phi_s)/(8*pi))*Iph1;
gu3_z=-((1/(4*pi))*Iz3 + (cos(2*phi_s)/(4*pi))*Iz1);

gu4_r=sin(phi_s)/(4*pi)*Ir2;
gu4_ph=-cos(phi_s)/(4*pi)*Iph2;
gu4_z=sin(phi_s)/(2*pi)*Iz2;
%%
cg1=-cos(lamda_s)*sin(delta_s);
cg2=cos(lamda_s)*cos(delta_s);
cg3=-sin(lamda_s)*sin(2*delta_s);
cg4=sin(lamda_s)*cos(2*delta_s);

u_r=cg1*gu1_r+cg2*gu2_r+cg3*gu3_r+cg4*gu4_r;
u_phi=cg1*gu1_ph+cg2*gu2_ph+cg3*gu3_ph+cg4*gu4_ph;
u_z=cg1*gu1_z+cg2*gu2_z+cg3*gu3_z+cg4*gu4_z;
t_vect=0:0.01:20;

T=0.4;
A=8/T^2;
[~,Domega]=fns_Source.get_Domega(A,T,omega_vect);
M0units=M_0/1e7/1e15;
% s_vect=pi*A*T*(1+exp(1i*omega_vect*T))./(pi^2-omega_vect.^2*T^2);
% idx_s_vect=find(abs(omega_vect-pi/T)<1e-12);
% s_vect(idx_s_vect)=1i*A*T/2;
% figure
% plot(omega_vect,s_vect)
Vphi_t=zeros(1,length(t_vect));
Vr_t=zeros(1,length(t_vect));
Vz_t=zeros(1,length(t_vect));
for i=1:length(t_vect)
    t=t_vect(1,i);
    %...................................................%
    Int_Ur=M0units*Domega.*u_r.*exp(-1i*omega_vect*t);
    Int_Vr=1i*omega_vect.*Int_Ur;
    Vr=exp(o_I*t)/pi*(trapz(omega_vect,Int_Vr));
    Vr_t(:,i)=Vr.';
    %...................................................%
    Int_Uphi=M0units*Domega.*u_phi.*exp(-1i*omega_vect*t);
    Int_Vphi=1i*omega_vect.*Int_Uphi;
    Vphi=exp(o_I*t)/pi*(trapz(omega_vect,Int_Vphi));
    Vphi_t(:,i)=Vphi.';
    %...................................................%
    Int_Uz=M0units*Domega.*u_z.*exp(-1i*omega_vect*t);
    Int_Vz=1i*omega_vect.*Int_Uz;
    Vz=exp(o_I*t)/pi*(trapz(omega_vect,Int_Vz));
    Vz_t(:,i)=Vz.';
    % integrand_hmg=s_vect.*Iz1.*exp(-1i*omega_vect*t)/2/pi;
    % I_hmg=exp(o_I*t)/pi*real(trapz(omega_vect,integrand_hmg));
    % u_zz_t_vect_hmg(:,i)=I_hmg.';
    display(strcat('t=',num2str(t),'s'));
end

%% Wald validation
factr=100;
figure;
subplot(3,1,1);
plot(t_vect,Vr_t*factr)
% ylim([-4 4])
subplot(3,1,2);
plot(t_vect,Vphi_t*factr)
% ylim([-5 5])
subplot(3,1,3);
plot(t_vect,Vz_t*factr)
% ylim([-1 1])
% figure;
% plot(t_vect,u_zz_t_vect_hmg)