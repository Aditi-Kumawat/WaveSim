clear;clc;close all

% Initial setup
r_vect=[0.5 1.5 3.3 5.6 7.5 10];
phi_s=deg2rad([70]);          %strike
delta_s=deg2rad([60]);        %dip
lamda_s=deg2rad([-20]); % a vector of slip or rake values
varied_angle=lamda_s;
M_L=2.1;
M_w=0.67 + 0.56*M_L + 0.046*M_L^2;
M_0=10^(1.5*(M_w+10.7));

% Frequency range
f_vect = 0.01:0.01:20;
omega_vect = 2 * pi * f_vect;
o_I = 1e-3;
t_vect=0:0.01:20;
nT=length(t_vect);
%% Source Parameters
T=0.1;
A=8/T^2;
S=1;
% [~,Domega]=fns_Source.get_Domega(A,T,omega_vect);
M0units=M_0/1e7/1e15;
Domega=pi*S*T*(1+exp(1i*omega_vect*T))./(pi^2-omega_vect.^2*T^2);
% idx_s_vect=find(abs(omega_vect-pi/T)<1e-12);
% Domega(idx_s_vect)=1i*S*T/2;
% figure
% plot(omega_vect,s_vect)
%% North south thingy
phi_az=30;
c_az=cos(phi_az);
s_az=sin(phi_az);

%%
%-----------------------------------%
soil_medium='Wald1996_LHS';
d_J=2.75;
dir=fns_inptMatPara1.form_dir(soil_medium,d_J);
%-----------------------------------%
Vi1_t_values = zeros(length(r_vect), nT);  % assuming t_vect length is 2001
Vi2_t_values = zeros(length(r_vect), nT);
Vz_t_values = zeros(length(r_vect), nT);
max_vt=zeros(3,length(r_vect));
for idx=1:length(r_vect)
    r=r_vect(idx);

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
    gu1_i1=gu1_r*c_az+gu1_ph*(-s_az);
    gu1_i2=gu1_r*s_az+gu1_ph*(c_az);

    gu2_i1=gu2_r*c_az+gu2_ph*(-s_az);
    gu2_i2=gu2_r*s_az+gu2_ph*(c_az);

    gu3_i1=gu3_r*c_az+gu3_ph*(-s_az);
    gu3_i2=gu3_r*s_az+gu3_ph*(c_az);

    gu4_i1=gu4_r*c_az+gu4_ph*(-s_az);
    gu4_i2=gu4_r*s_az+gu4_ph*(c_az);


    %%
    cg1=-cos(lamda_s)*sin(delta_s);
    cg2=cos(lamda_s)*cos(delta_s);
    cg3=-sin(lamda_s)*sin(2*delta_s);
    cg4=sin(lamda_s)*cos(2*delta_s);

    u_i1=cg1*gu1_i1+cg2*gu2_i1+cg3*gu3_i1+cg4*gu4_i1;
    u_i2=cg1*gu1_i2+cg2*gu2_i2+cg3*gu3_i2+cg4*gu4_i2;
    u_z=cg1*gu1_z+cg2*gu2_z+cg3*gu3_z+cg4*gu4_z;

    Vi2_t=zeros(1,nT);
    Vi1_t=zeros(1,nT);
    Vz_t=zeros(1,nT);
    for i=1:nT
        t=t_vect(1,i);
        %...................................................%
        Int_Ui1=M0units*Domega.*u_i1.*exp(-1i*omega_vect*t);
        Int_Vi1=1i*omega_vect.*Int_Ui1;
        Vi1=exp(o_I*t)/pi*(trapz(omega_vect,Int_Vi1));
        Vi1_t(:,i)=Vi1.';
        
        %.............0.1 ......................................%
        Int_Ui2=M0units*Domega.*u_i2.*exp(-1i*omega_vect*t);
        Int_Vi2=1i*omega_vect.*Int_Ui2;
        Vi2=exp(o_I*t)/pi*(trapz(omega_vect,Int_Vi2));
        Vi2_t(:,i)=Vi2.';
        %...................................................%
        Int_Uz=M0units*Domega.*u_z.*exp(-1i*omega_vect*t);
        Int_Vz=1i*omega_vect.*Int_Uz;
        Vz=exp(o_I*t)/pi*(trapz(omega_vect,Int_Vz));
        Vz_t(:,i)=Vz.';
    end
    Vi1_t_values(idx, :) = Vi1_t;
    Vi2_t_values(idx, :) = Vi2_t;
    Vz_t_values(idx, :) = Vz_t;
    max_vt(1,idx)=abs(max(Vi1_t));
    max_vt(2,idx)=abs(max(Vi2_t));
    max_vt(3,idx)=abs(max(Vz_t));
end

% Plotting results
factr=1e3;
y_lbl='HHE (mm/s)';
filename = ['Poing_HHE_ML_', num2str(M_L), '.pdf'];
fns_PlotSimulatedFreefield.plot_HH_ENZ(t_vect,Vi1_t_values,factr,r_vect,y_lbl,filename);

y_lbl='HHN (mm/s)';
filename = ['Poing_HHN_ML_', num2str(M_L), '.pdf'];
fns_PlotSimulatedFreefield.plot_HH_ENZ(t_vect,Vi2_t_values,factr,r_vect,y_lbl,filename);

y_lbl='HHZ (mm/s)';
filename = ['Poing_HHZ_ML_', num2str(M_L), '.pdf'];
fns_PlotSimulatedFreefield.plot_HH_ENZ(t_vect,Vz_t_values,factr,r_vect,y_lbl,filename);

fns_PlotSimulatedFreefield.plot_maxVwithr(factr,r_vect,max_vt);