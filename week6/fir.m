clc;
clear;
close all;

%% === User Inputs ===
fprintf('--- FIR Filter Design using Window Method ---\n');

passripple = input('Enter passband ripple: ');
stopatt = input('Enter stopband attenuation: ');
wp = input('Enter passband frequency: ');
ws = input('Enter stopband frequency: ');

%% === Estimate Filter Order ===
deltaw = abs(ws - wp);
if deltaw <= 0
    error('Stopband frequency must be greater than passband frequency.');
end

Nest = ceil((-10 * log10(passripple * stopatt) - 15) / (14 * deltaw));
if mod(Nest, 2) == 0
    Nest = Nest + 1;
end

fprintf('Estimated Filter Order: N = %d\n', Nest);

tau = (Nest - 1) / 2;
n = 0:Nest - 1;

%% === Filter and Window Types ===
filttypes = {'lowpass', 'highpass', 'bandpass', 'bandstop'};
wintypes = {'rectangular', 'modifiedrect', 'bartlett', 'hanning', 'hamming', 'blackman'};

colors = lines(length(wintypes));

for f = 1:length(filttypes)
    num_windows = length(wintypes);
    hdcell = cell(1, num_windows);
    freqcell = cell(1, num_windows);
    wincell = cell(1, num_windows);

    % Precompute impulse & frequency responses
    for widx = 1:num_windows

        % --- Ideal Impulse Response ---
        switch filttypes{f}
            case 'lowpass'
                hd = hdlpf(n, tau, wp);
            case 'highpass'
                hd = hdhpf(n, tau, wp);
            case 'bandpass'
                if wp >= ws
                    error('Lower cutoff wp must be < upper cutoff ws.');
                end
                hd = hdbpf(n, tau, wp, ws);
            case 'bandstop'
                if wp >= ws
                    error('Lower cutoff wp must be < upper cutoff ws.');
                end
                hd = hdbsf(n, tau, wp, ws);
        end

        % --- Window selection ---
        switch wintypes{widx}
            case 'rectangular'
                win = rectwin(n, Nest);
            case 'modifiedrect'
                win = modrectwin(n, Nest);
            case 'bartlett'
                win = bartlettwin(n, Nest);
            case 'hanning'
                win = hanningwin(n, Nest);
            case 'hamming'
                win = hammingwin(n, Nest);
            case 'blackman'
                win = blackmanwin(n, Nest);
        end

        hfir = hd .* win;

        % --- Frequency Response (manual FFT) ---
        Hf = manualfft(hfir, 2048);
        freqaxis = linspace(-pi, pi, length(Hf)) / pi * 2; % from -2 to 2

        freqcell{widx} = Hf;
        hdcell{widx} = hfir;
        wincell{widx} = win;
    end

    %% --- Plotting (6x3 layout: Phase, Magnitude, Window) ---
    figure('Name', [filttypes{f}, ' Filter Responses'], 'NumberTitle', 'off');
    sgtitle([filttypes{f}, ' Filter Responses']);

    for widx = 1:num_windows
        phase_idx = (widx - 1) * 3 + 1;
        mag_idx = (widx - 1) * 3 + 2;
        win_idx = (widx - 1) * 3 + 3;

        freqaxis = linspace(-pi, pi, length(freqcell{widx})) / pi * 2;

        % === Phase Response ===
        subplot(6, 3, phase_idx);
        plot(freqaxis, unwrapmanualphase(freqcell{widx}), ...
            'Color', colors(widx, :), 'LineWidth', 1.5);
        title([wintypes{widx}, ' Phase'], 'Interpreter', 'none');
        xlabel('\omega / \pi');
        ylabel('Phase (rad)');
        grid on;
        xlim([-2 2]);

        % === Magnitude Response ===
        subplot(6, 3, mag_idx);
        plot(freqaxis, manualabs(freqcell{widx}), ...
            'Color', colors(widx, :), 'LineWidth', 2);
        title(sprintf('%s |H(\\omega)|', wintypes{widx}), 'Interpreter', 'none');
        xlabel('\omega / \pi');
        ylabel('Magnitude');
        grid on;
        ylim([0 1.2 * max(manualabs(freqcell{widx}))]);
        xlim([-2 2]);

        % === Window Function ===
        subplot(6, 3, win_idx);
        stem(n, wincell{widx}, 'filled', 'Color', colors(widx, :));
        title([wintypes{widx}, ' Window'], 'Interpreter', 'none');
        xlabel('n');
        ylabel('w[n]');
        grid on;
    end
end

%% ======= Ideal Impulse Response Functions =======
function hd = hdlpf(n, tau, wp)
hd = zeros(size(n));
for i = 1:length(n)
    if n(i) ~= tau
        hd(i) = sin(wp * (n(i) - tau)) / (pi * (n(i) - tau));
    else
        hd(i) = wp / pi;
    end
end
end

function hd = hdhpf(n, tau, wp)
hd = zeros(size(n));
for i = 1:length(n)
    if n(i) ~= tau
        hd(i) = -sin(wp * (n(i) - tau)) / (pi * (n(i) - tau));
    else
        hd(i) = 1 - wp / pi;
    end
end
end

function hd = hdbpf(n, tau, wp1, wp2)
hd = zeros(size(n));
for i = 1:length(n)
    if n(i) ~= tau
        hd(i) = (sin(wp2 * (n(i) - tau)) - sin(wp1 * (n(i) - tau))) / (pi * (n(i) - tau));
    else
        hd(i) = (wp2 - wp1) / pi;
    end
end
end

function hd = hdbsf(n, tau, wp1, wp2)
hd = zeros(size(n));
for i = 1:length(n)
    if n(i) ~= tau
        hd(i) = (sin(wp1 * (n(i) - tau)) - sin(wp2 * (n(i) - tau))) / (pi * (n(i) - tau));
    else
        hd(i) = 1 - (wp2 - wp1) / pi;
    end
end
end

%% ======= Window Functions =======
function w = rectwin(n, N)
w = ones(size(n));
end

function w = modrectwin(n, N)
w = ones(size(n));
w(1) = 0.5;
w(end) = 0.5;
end

function w = bartlettwin(n, N)
w = zeros(size(n));
for i = 1:length(n)
    if n(i) <= (N - 1) / 2
        w(i) = 2 * n(i) / (N - 1);
    else
        w(i) = 2 - 2 * n(i) / (N - 1);
    end
end
end

function w = hanningwin(n, N)
w = 0.5 * (1 - cos(2 * pi * n / (N - 1)));
end

function w = hammingwin(n, N)
alpha = 0.54;
w = alpha - (1 - alpha) * cos(2 * pi * n / (N - 1));
end

function w = blackmanwin(n, N)
w = 0.42 - 0.5 * cos(2 * pi * n / (N - 1)) + 0.08 * cos(4 * pi * n / (N - 1));
end

%% ======= Manual FFT =======
function X = manualfft(x, N)
x = [x, zeros(1, N - length(x))];
X = zeros(1, N);
for k = 0:N - 1
    sumval = 0;
    for n = 0:N - 1
        sumval = sumval + x(n + 1) * exp(-1j * 2 * pi * k * n / N);
    end
    X(k + 1) = sumval;
end
X = fftshiftmanual(X);
end

function Xs = fftshiftmanual(X)
N = length(X);
p = floor(N / 2);
Xs = [X(p + 1:end), X(1:p)];
end

%% ======= Manual Magnitude and Phase =======
function mag = manualabs(X)
mag = sqrt(real(X).^2 + imag(X).^2);
end

function ph = unwrapmanualphase(X)
ph = anglemanual(X);
for i = 2:length(ph)
    while ph(i) - ph(i - 1) > pi
        ph(i:end) = ph(i:end) - 2 * pi;
    end
    while ph(i) - ph(i - 1) < -pi
        ph(i:end) = ph(i:end) + 2 * pi;
    end
end
end

function ph = anglemanual(X)
ph = atan2(imag(X), real(X));
end
