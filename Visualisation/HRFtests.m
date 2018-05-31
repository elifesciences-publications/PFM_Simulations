%Looks at some of the properties of the HRF
% Sam Harrison

clear all; close all; clc
addpath(genpath('~/Documents/Code/MATLAB/'))
prettyFigures();

%% 

dt = 1e-3; %Sampling rate
T = 30;     %Total time

t = linspace(0, T, 1+T/dt);

%% Load the HRF

HRF = load('/usr/local/fsl/etc/default_flobs.flobs/hrfbasisfns.txt');
HRF = HRF(:,1);
dtHRF = 0.05;   %FLOBS sampling rate

%Interpolate HRF to requested rate
tHRF = dtHRF*(1:length(HRF)) - dtHRF;
tNew = dt * (0:floor(tHRF(end)/dt));
HRF = interp1(tHRF, HRF, tNew, 'cubic', 0);

%Set to DC gain of 1
HRF = HRF / trapz(tNew, HRF);
%HRF = HRF / sum(HRF);

%Plot
figure; plot(tNew, HRF); xlim([tNew(1) tNew(end)])
xlabel('Time (s)'); title('HRF')

%% Do some tests!

%Block
x = zeros(length(t),1);
x( (t>=2.5) & (t<=7.5) ) = 1;

y = dt * conv(x, HRF);

figure; plot(t, x); hold on; plot(dt*(1:length(y))-dt, y, 'r')
ylim([-0.2 1.2]); xlabel('Time (s)')

legend('Signal', 'HRF(signal)')


%Block wave
x = zeros(length(t),1);
x( (t>=2.5) & (t<=7.5) ) = 1;
x( (t>=2.5) & (t<=7.5) ) = sin( 2*pi*10*(dt*(1:sum(x))-dt) );

y = dt * conv(x, HRF);
z = dt * conv(abs(x), HRF);

figure; plot(t, x); hold on; 
plot(dt*(1:length(y))-dt, y, 'g')
plot(dt*(1:length(z))-dt, z, 'r')
ylim([-1.2 1.2]); xlabel('Time (s)')

legend('Signal', 'HRF(signal)', 'HRF(abs(signal))')



%Envelope wave
x = sin( 2*pi*10*t ) .* sin( 2*pi*0.063*t );

y = dt * conv(x, HRF);
z = dt * conv(abs(x), HRF);

figure; plot(t, x); hold on; 
plot(dt*(1:length(y))-dt, y, 'g')
plot(dt*(1:length(z))-dt, z, 'r')
ylim([-1.2 1.2]); xlabel('Time (s)')

legend('Signal', 'HRF(signal)', 'HRF(abs(signal))')

figure; plot(t, x); hold on; 
plot(dt*(1:length(y))-dt, y, 'g')
plot(dt*(1:length(z))-dt, z, 'r')
ylim([-1.2 1.2]); xlim([12 17]); xlabel('Time (s)')

legend('Signal', 'HRF(signal)', 'HRF(abs(signal))')
