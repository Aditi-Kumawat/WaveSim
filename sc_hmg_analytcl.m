% Pay attention to k_vect range while changing the frequency range
%% Medium properties
clear;clc
rho= 2.3;
alpha=2.4; % P-wave velocity in m/s
beta=1; % S-wave velocity in m/s
mu = rho .* beta.^2;
% rho=mu/beta^2; %medium density kg/m^3

%% K_mat chunks
% Division of k_vect: to make the formulation faster
dk = 1.25e-4;
kmin=1.25e-4;
kmax=20;
N_laps =2000;
k_vect=0:dk:kmax;
k_N_1=round(length(k_vect)/N_laps);
k_N_vect=zeros(1,N_laps);
k_N_vect(1,1)=1;
for i_k_N=1:N_laps-1
    k_N_vect(1,i_k_N+1)=i_k_N*k_N_1+1;
end
k_N_vect=[k_N_vect length(k_vect)];

k_mat=zeros(N_laps,k_N_1+2);
for i_laps=1:N_laps
    k_1=k_vect(k_N_vect(i_laps));
    k_2=k_vect(k_N_vect(i_laps+1));
    k_vect_1=k_1:dk:k_2;
    k_mat(i_laps,1:length(k_vect_1))=k_vect_1;
end
k_mat_final=k_mat(1:end-1,1:end-1); % the last one entry is zero
k_vect_last=k_mat(end,:);
idx_last_0=find(k_vect_last<1e-12);
k_vect_last(idx_last_0)=[];

%% f_hz_range
% f_hz_vect=5;
f_hz_vect=0.01:0.01:5;
omega_vect=2*pi*f_hz_vect;
o_I=1e-3;
tic;
x=10; % (km)
%% Evaluation of frequency domain green's function
g_z_0_omega_mat=zeros(N_laps-1,length(omega_vect));

for i_laps=1:N_laps-1
    display(strcat('lap=',num2str(i_laps)));
    k_vect_1=k_mat_final(i_laps,:);
    J0= besselj(0,k_vect_1*x); % Bessel function
    for i_omg_1=1:length(omega_vect)
        omega=omega_vect(i_omg_1); 
        %         display(omega)
        k_beta_vect_1=(omega+1i*o_I)./beta;
        k_alpha_vect_1=(omega+1i*o_I)./alpha;
        
        gamma_vect_1=(k_vect_1.^2-k_beta_vect_1.^2).^(1/2);
        nu_vect_1=(k_vect_1.^2-k_alpha_vect_1.^2).^(1/2);
        
        chi_vect_1=k_vect_1.^2+gamma_vect_1.^2;
        nu_gamma_vect_1=(gamma_vect_1.*nu_vect_1);
        
        n_f_RR_vect_1=nu_vect_1.*(chi_vect_1-2*k_vect_1.^2);
        Denom_vect_1=mu*(4*k_vect_1.^2.*nu_gamma_vect_1-chi_vect_1.^2);
        
        f_RR_vect_1=n_f_RR_vect_1./Denom_vect_1; % f_RR
        
        g_z_int_0_vect_1=k_vect_1.*f_RR_vect_1.*J0;
        g_z_0_omega=(trapz(k_vect_1,g_z_int_0_vect_1));
        g_z_0_omega_mat(i_laps,i_omg_1)=g_z_0_omega;
    end
end
g_z_0_omega_vect_1=sum(g_z_0_omega_mat);
g_z_0_omega_mat_idx=isnan(g_z_0_omega_mat);
g_z_0_omega_mat_nan=g_z_0_omega_mat(g_z_0_omega_mat_idx);

% for the last k_vect
g_z_0_omega_vect_last=zeros(1,length(omega_vect));
k_x_vect_last=(k_vect_last)*x;
J0_last=besselj(0,k_vect_last*x);
display(strcat('lap=',num2str(N_laps)));

for i_omg_2=1:length(omega_vect)
    omega=omega_vect(i_omg_2);
    k_beta_vect_last=(omega+1i*o_I)./beta;
    k_alpha_vect_last=(omega+1i*o_I)./alpha;
    
    gamma_vect_last=(k_vect_last.^2-k_beta_vect_last.^2).^(1/2);
    nu_vect_last=(k_vect_last.^2-k_alpha_vect_last.^2).^(1/2);
    
    chi_vect_last=k_vect_last.^2+gamma_vect_last.^2;
    nu_gamma_vect_last=(gamma_vect_last.*nu_vect_last);
    
    n_f_RR_vect_last=nu_vect_last.*(chi_vect_last-2*k_vect_last.^2);
    Denom_vect_last=mu*(4*k_vect_last.^2.*nu_gamma_vect_last-chi_vect_last.^2);
    
    f_RR_vect_last=n_f_RR_vect_last./Denom_vect_last;
    
    g_z_int_0_vect_last=k_vect_last.*f_RR_vect_last.*J0_last;
    g_z_0_omega_vect_last(i_omg_2)=(trapz(k_vect_last,g_z_int_0_vect_last));
end

%% frequency domain green's function
g_z_0_omega_vect=g_z_0_omega_vect_1+g_z_0_omega_vect_last;
figure
plot(f_hz_vect,g_z_0_omega_vect)
%% half sine pulse
t_vect=0:1e-3:20;
T=1;
A=1;
s_vect=pi*A*T*(1+exp(1i*omega_vect*T))./(pi^2-omega_vect.^2*T^2);
idx_1=find(abs(omega_vect-pi/T)<1e-12);
s_vect(idx_1)=1i*A*T/2;

%% Inverse Fourier Transform
int_g_z_0_hsp_vect=g_z_0_omega_vect.*s_vect;
g_z_0_t_vect=zeros(1,length(t_vect));
for i_t=1:length(t_vect)
    t=t_vect(1,i_t);
    integrand=int_g_z_0_hsp_vect.*exp(-1i*omega_vect*t);
    g_z_0_t_vect(1,i_t)=trapz(omega_vect,integrand,2);
    display(strcat('t=',num2str(t),'s'));
    
end
% figure
% plot(omega_vect,int_g_z_0_hsp_vect)
% int_g_z_0_hsp_mat=(ones(length(t_vect),1)*int_g_z_0_hsp_vect).*exp_mat;
% g_z_0_t=trapz(omega_vect,int_g_z_0_hsp_mat,2);
e_o_I_vect=exp(o_I*t_vect)*(1/pi);
g_z_0_t_vect=(1/2/pi)*g_z_0_t_vect.*e_o_I_vect; % 1/pi for fourier transform and 1/pi/mu for inverse transform from wavenumber domain

t_ex=toc;
%%
ha_col=@colors;
figure
plot(t_vect,-g_z_0_t_vect,'color',ha_col('ball blue'),'LineStyle','-','LineWidth',1.5);
% ylim([-1.1 0.1])
% xlim([0.2 2])
xlabel({'time~(s)';'\bf{(a)}'},'FontSize',11,'Interpreter','latex')
ylabel('Normalized~Deflection,~${w}_0({\hat{z}})$','FontSize',11,'Interpreter','latex')
% legend({strcat('$x=1m$')},...
%     'Location','southeast','NumColumns',4,...
%     'FontSize',10,'Interpreter','latex','box','off')
ax.XTickLabelMode = 'auto';
ax.YTickLabelMode = 'auto';

set(gcf,'Units','inches', 'Position', [4 3 5 2.5]);
grid off
set(gcf,'renderer','Painters')