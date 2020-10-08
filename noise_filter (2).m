clear;
load('09_01_dt.mat');
load('09_01_signals.mat');
load('11_01_dt.mat');
load('11_01_signals.mat');

%% constants and variables
% parameters
DT = 0.005; % time between measurements
x  =  x12_11_01;  % raw data signal - vector

% constants
RD2V = 5e3 / 1024; % convertion factor from Arduino raw data to mV
OFFSET = 350;  % Arduino offset in mV (?)
A = 4.8;  % chip Amplification - (estimated  for now)

N = length(x);

% time
t_min = 0; % just for convenience
t_max = t_min + (N*DT);
t = linspace(t_min, t_max, N); % time vector

% frequencies
Fs = 1 / DT;
f = linspace(0, Fs / 2, N/2 + 1);  % frequency vector
% omega = f .* (2*pi);
samples_per_freq = N/Fs;

% plot raw data
figure;
subplot(2, 2, 1)
plot(t, x)
grid on
title("Raw Data Signal")
xlabel("t (sec)")
ylabel("V (Arduino)")

% basic processing
x = x .* (RD2V/A) - OFFSET; % signal without amplification

% After basic processing
subplot(2, 2, 2)
plot(t, x)
grid on
title("original signal")
xlabel("t (sec)")
ylabel("V (mV)")


%% frequency spectrum

X = abs(fftshift(fft(x)));
X = X(end-(N/2):end);

subplot(2, 2, 3)
plot(f, X);
grid on
title("signal spectrum")
xlabel("f (Hz)")
ylabel("V (mV/Hz)")

%% original Vs filtered signal 

% Band-pass filtering for frog heart using:
%       Fstop1 = 0.2, Fpass1 = 0.4, Fpass2 = 5, Fstop2 = 6
B = [0.028  0.053 0.071  0.053 0.028];  % Numerator
A = [1.000 -2.026 2.148 -1.159 0.279];  % Denominator
y = filter(B, A, x);

% After
X = fftshift(fft(x));
cut_off = 8;  % cut_off frequency in Hz
wind_size = round(cut_off*samples_per_freq);  % wind size just for positive frequncies
filter = zeros(length(X), 1);

% wind = rectwin(2*wind_size+1);  % Naive window
wind = kaiser(2*wind_size+1);

filter((end/2)-wind_size:(end/2)+wind_size) = wind;
X = X .* filter;

subplot(2, 2, 4)
plot(t, ifft(ifftshift(X)));
% plot(t, y);
grid on
title("filtered signal")
xlabel("t (sec)")
ylabel("V (mV)")


