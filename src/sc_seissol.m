clear;close all;clc;
% Import the data from the 'seissol' folder
data = readtable('seissol/XX0303.csv');

% Extract columns based on provided header names
time = data.Time_s_;
z_velocity = data.HHZ_m_s_;
north_velocity = data.HHN_m_s_;
east_velocity = data.HHE_m_s_;

% Plotting
figure;
subplot(3,1,1);
plot(time, z_velocity, 'b', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('HHZ [m/s]');
xlim([0,5.5])

subplot(3,1,2);
plot(time, north_velocity, 'r', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('HHN [m/s]');
xlim([0,5.5])

subplot(3,1,3);
plot(time, east_velocity, 'g', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('HHE [m/s]');
xlim([0,5.5])

% Given coordinates
epiC = [4466110, 5331600];
XX0303 = [4464460, 5329950];

% Calculate the angle in degrees
angle_deg = atan2d(XX0303(2) - epiC(2), XX0303(1) - epiC(1));

% Calculate the distance between the two points
distance = sqrt((XX0303(1) - epiC(1))^2 + (XX0303(2) - epiC(2))^2);

% Plotting
figure;
plot([epiC(1), XX0303(1)], [epiC(2), XX0303(2)], '-o');
text(epiC(1), epiC(2), 'Epicenter', 'VerticalAlignment', 'bottom');
text(XX0303(1), XX0303(2), 'Receiver', 'VerticalAlignment', 'bottom');
midpoint = (epiC + XX0303) / 2;
text(midpoint(1), midpoint(2), ['Distance: ', num2str(distance, '%.2f')], 'VerticalAlignment', 'top');
title(['Angle = ', num2str(angle_deg), '°']);
xlabel('X Coordinate');
ylabel('Y Coordinate');
grid on;
axis equal;

