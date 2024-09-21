% Audio Signal Processing Script
% This MATLAB script performs various audio signal processing tasks including 
% time-domain and frequency-domain analysis, signal segment extraction, 
% cross-correlation, resampling, Fourier Transform, custom DFT calculation, 
% noise addition, modulation, demodulation, and filtering.
%
% Author: Ali Pazani

% Clearing the workspace
clc; clear; close all;

% Load the audio file
[x, fs] = audioread('ghayegh.wav');
x = x(:,1); % Use the first channel if stereo
N = length(x);
t = (0:N-1) / fs;

%% Part 1: Plot the Original Signal in the Time Domain
figure;
plot(t, x);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Original "ghayegh" Signal in Time Domain');

%% Part 2: Extract and Plot a Specific Segment of the Signal
a = 63; % Start time in seconds
segment_start = a * fs;
segment_end = (a + 10) * fs;
y = x(segment_start:segment_end); % Extract 10-second segment
n = (0:length(y)-1) / fs + a; % Time vector for the segment

figure;
plot(n, y);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Extracted "ghayegh2" Signal in Time Domain');
axis([a, a + 10, -1, 1]); % Axis limits
audiowrite('ghayegh2.wav', y, fs);

%% Part 3: Cross-Correlation of Two Signal Segments
y1 = y(1:0.2*fs); % First 200ms segment
y2 = y(0.1*fs:0.2*fs); % Segment from 100ms to 200ms
[r, lags] = xcorr(y1, y2); % Cross-correlation
timeLags = lags / fs * 1000; % Convert to milliseconds

figure;
subplot(3,1,1);
plot(linspace(0, 20, length(y1)), y1);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Segment y[0:20ms]');

subplot(3,1,2);
plot(linspace(10, 20, length(y2)), y2);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Segment y[10:20ms]');

subplot(3,1,3);
plot(timeLags, r);
xlabel('Time Lag (ms)');
ylabel('Amplitude');
title('Cross-Correlation of y[0:20ms] and y[10:20ms]');

%% Part 4: Resampling the Signal (Downsampling and Upsampling)
y_downsampled = y(1:2:end); % Downsample by a factor of 2 (y[2n])
y_upsampled = zeros(2*length(y), 1); % Upsample by inserting zeros (y[n/2])
y_upsampled(1:2:end) = y;

% Time vectors for plotting
t_downsampled = linspace(0, 10, length(y_downsampled));
t_upsampled = linspace(0, 10, length(y_upsampled));

figure;
subplot(3,1,1);
plot(linspace(0, 10, length(y)), y);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original "ghayegh2" Signal');

subplot(3,1,2);
plot(t_downsampled, y_downsampled);
xlabel('Time (s)');
ylabel('Amplitude');
title('Downsampled Signal y[2n]');

subplot(3,1,3);
plot(t_upsampled, y_upsampled);
xlabel('Time (s)');
ylabel('Amplitude');
title('Upsampled Signal y[n/2]');

% Save the modified signals
audiowrite('ghayegh2_1.wav', y_downsampled, fs);
audiowrite('ghayegh1_2.wav', y_upsampled, fs);

%% Part 5: Fourier Transform and Frequency Domain Representation
[y, fs] = audioread('ghayegh2.wav'); % Load the segment 'ghayegh2.wav'
N = length(y);
f = linspace(-pi, pi, N); % Frequency vector

Y = fftshift(fft(y)); % Compute and shift FFT
magnitude_Y = abs(Y);

figure;
plot(f, magnitude_Y);
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Magnitude Spectrum of "ghayegh2"');
grid on;

%% Part 6: Analyze Frequency Changes Due to Resampling
[y_downsampled, ~] = audioread('ghayegh2_1.wav'); % Load downsampled signal
[y_upsampled, ~] = audioread('ghayegh1_2.wav'); % Load upsampled signal

% Fourier Transform for the downsampled signal
N_downsampled = length(y_downsampled);
Y_downsampled = fftshift(fft(y_downsampled));
magnitude_Y_downsampled = abs(Y_downsampled);
f_downsampled = linspace(-pi, pi, N_downsampled);

% Fourier Transform for the upsampled signal
N_upsampled = length(y_upsampled);
Y_upsampled = fftshift(fft(y_upsampled));
magnitude_Y_upsampled = abs(Y_upsampled);
f_upsampled = linspace(-pi, pi, N_upsampled);

figure;
subplot(2,1,1);
plot(f_downsampled, magnitude_Y_downsampled);
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Magnitude Spectrum of Downsampled Signal y[2n] (ghayegh2_1)');
grid on;

subplot(2,1,2);
plot(f_upsampled, magnitude_Y_upsampled);
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Magnitude Spectrum of Upsampled Signal y[n/2] (ghayegh1_2)');
grid on;

%% Part 7: Compression and Frequency Analysis
downsample_factor = 3; % Downsample the signal by a factor of 3
y_downsampled = downsample(y, downsample_factor);
audiowrite('ghayegh7_3.wav', y_downsampled, fs/downsample_factor); % Save the downsampled signal

N_downsampled = length(y_downsampled);
f_downsampled = linspace(-pi, pi, N_downsampled);
Y_downsampled = fftshift(fft(y_downsampled));
magnitude_Y_downsampled = abs(Y_downsampled);

figure;
plot(f_downsampled, magnitude_Y_downsampled);
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Magnitude Spectrum of Downsampled "ghayegh7_3"');
grid on;

%% Part 8: Custom DFT Calculation
N_dft = 44;
DFT = zeros(N_dft, 1);
for k = 1:N_dft
    for n = 1:length(y)
        DFT(k) = DFT(k) + y(n) * exp(-1j * 2 * pi * (k-1) * (n-1) / N_dft);
    end
end
DFT = fftshift(DFT); % Center zero frequency component
magnitude_DFT = abs(DFT);

figure;
plot(1:N_dft, magnitude_DFT);
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Magnitude Spectrum using Custom DFT');
grid on;

%% Part 9: Adding White Gaussian Noise (AWGN) to the Signal
[x, fs] = audioread('ghayegh2.wav');
Time = (1/fs) * length(x);
T = linspace(0, Time, length(x));

out = awgn(x, 15); % Add AWGN (white Gaussian noise) with SNR = 15dB

figure;
plot(T, out);
hold on;
plot(T, x);
legend({'Signal with AWGN', 'Original Signal'}, 'Location', 'northwest');
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs Noisy "ghayegh-n" Signal');
axis([0, 10, -2, 2]);
audiowrite('ghayegh-n.wav', out, fs);

%% Part 10: Signal Modulation, Demodulation, and Filtering
% Modulation with a cosine signal
t = linspace(1, 10, length(x));
g = x' .* cos(2*pi*6500.*t);

audiowrite('y_modulated.wav', g, fs);

% Fourier Transform of the modulated signal
DFg = fft(g);
w = linspace(0, 2*pi, length(g));

figure;
plot(w, abs(DFg));
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Modulated Signal "y-modulated"');
axis([0, 2*pi, -Inf, Inf]);

% Demodulation
h = g .* cos(2*pi*6500.*t);
audiowrite('y_demodulated.wav', h, fs);

DFh = fft(h);

figure;
plot(w, abs(DFh));
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Demodulated Signal "y-demodulated"');
axis([0, 2*pi, -Inf, Inf]);

% Filtering the demodulated signal
load('filt_taps');
u = conv(h, filt_taps, 'same'); % Convolution of h with filt_taps
audiowrite('y-demodulated2.wav', u, fs);

DFu = fft(u);

figure;
plot(w, abs(DFu));
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
title('Filtered Demodulated Signal "Convolution Result"');
axis([0, 2*pi, -Inf, Inf]);

