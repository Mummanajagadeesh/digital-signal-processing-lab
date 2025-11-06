clc;
clear all;
close all;

% ===== Continuous-time sine wave =====
t = 0:0.0001:0.1;
A = 10;
f = 50;

subplot(3, 3, 1);
x = A * sin(2 * pi * f * t);
plot(t, x, 'LineWidth', 1.2);
ylim([-20 20]);
grid on;
title('Continuous-Time Sine Wave (50 Hz)');
xlabel('Time (s)');
ylabel('Amplitude');

% ===== Sampling at different frequencies =====
fs_values = [50 100 200 400 500 1000 2000 4000];

for i = 1:length(fs_values)
    fs = fs_values(i);
    ts = 1 / fs;
    n = 0:ts:0.1;
    x_samp = A * sin(2 * pi * f * n);

    subplot(3, 3, i + 1);
    stem(n, x_samp, 'filled');
    ylim([-20 20]);
    grid on;
    title(['Sampled at f_s = ' num2str(fs) ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
end
