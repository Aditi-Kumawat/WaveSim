clear;close all;clc
% Parameters
T = 0.4;          % Total period of the waveform
A0 = 8/T^2;           % Amplitude of the waveform
% t = linspace(0, T, 1000);  % Time vector from 0 to T
dt=1e-3;
t = 0:dt:T;

% Compute the waveform
waveform = arrayfun(@(x) fns_Source.double_triangle_wave(x, T, A0), t);

% Numerical integration of the waveform to get velocity
velocity = cumtrapz(t, waveform);

% Numerical integration of the velocity to get displacement
displacement = cumtrapz(t, velocity);

%%
% Plot the waveform (acceleration), velocity, and displacement
t_ex=T+dt:dt:10;
fn_ex=zeros(1,length(t_ex));
t=[t t_ex];
waveform=[waveform fn_ex];
velocity=[velocity fn_ex];
displacement=[displacement fn_ex];
figure;

subplot(3, 1, 1);
plot(t, waveform);
xlabel('Time');
ylabel('Acceleration');
title('Double Triangular Wave - Acceleration');
grid on;

subplot(3, 1, 2);
plot(t, velocity);
xlabel('Time');
ylabel('Velocity');
title('Velocity');
grid on;

subplot(3, 1, 3);
plot(t, displacement);
xlabel('Time');
ylabel('Displacement');
title('Displacement');
grid on;
%%
fn_in= waveform;
Fs=1/dt;
nfft = 2^nextpow2(length(t));
freq = Fs / 2 * linspace(0, 1, nfft/2+1);
fn_fft = fft(fn_in, nfft) * (1/Fs);
%             fn_fft_ss = 2 * fn_fft(1:nfft/2+1,:); %% orignal
fn_fft_ss =fn_fft(1:nfft/2+1);

%%
[Aomega,Domega]=fns_Source.get_Domega(A0,T,freq.*2*pi);
%%
figure
plot(freq,fn_fft_ss)
hold on
plot(freq,Aomega,'-.')
xlim([0 50])
figure
plot(freq,Domega,'k-.')
xlim([0 50])
%%

