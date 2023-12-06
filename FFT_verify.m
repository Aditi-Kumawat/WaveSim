close all;clear;clc
f_hz_vect=0.01:0.01:20;
omega_vect=2*pi*f_hz_vect;
T=0.1;
A=1;
s_vect=pi*A*T*(1+exp(1i*omega_vect*T))./(pi^2-omega_vect.^2*T^2);
idx_1=find(abs(omega_vect-pi/T)<1e-12);
s_vect(idx_1)=1i*A*T/2;

%%
t_vect=-0.5:1e-3:50;
rect_vect=rectangularPulse(0, T, t_vect);
cos_vect=A*sin(pi/T*t_vect);
g_t_vect=cos_vect.*rect_vect;
figure
plot(t_vect,g_t_vect)
%%
dt=1e-3;
fft_unscaled=fft(g_t_vect);

sample_rate=1/dt;
Nsamples=length(t_vect);
df=sample_rate./Nsamples;
freq1=(0:df:sample_rate-df)';
figure
plot(freq1,fft_unscaled.*(1./sample_rate))
hold on
plot(f_hz_vect,s_vect,'-.')
xlim([0,20])
%%
%%
fn_in= g_t_vect;
Fs=1/dt;
nfft = 2^nextpow2(length(t_vect));
freq = Fs / 2 * linspace(0, 1, nfft/2+1);
fn_fft = fft(fn_in, nfft) * (1/Fs);
%             fn_fft_ss = 2 * fn_fft(1:nfft/2+1,:); %% orignal
fn_fft_ss =fn_fft(1:nfft/2+1);


figure
plot(f_hz_vect,s_vect)
hold on
plot(freq1,fft_unscaled.*(1./sample_rate),'-.')
hold on
plot(freq,fn_fft_ss,'--')
xlim([0 20])