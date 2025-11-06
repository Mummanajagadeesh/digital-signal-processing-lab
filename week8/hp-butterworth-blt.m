clc;
clear;

omega_p = input('Enter passband edge frequency (rad/s): ');
omega_s = input('Enter stopband edge frequency (rad/s): ');
Ap = input('Enter passband ripple (in dB): ');
As = input('Enter stopband attenuation (in dB): ');
fs = input('Enter sampling frequency (in Hz): ');

T = 1 / fs;

omega_p = (2 / T) * tan(omega_p * T / 2);
omega_s = (2 / T) * tan(omega_s * T / 2);

temp = omega_p;
omega_p = omega_s;
omega_s = temp;

del_p = 10^(-Ap / 20);
del_s = 10^(-As / 20);

a = (1 / del_s^2) - 1;
b = (1 / del_p^2) - 1;
c = log10(omega_s / omega_p);
x = log10(sqrt(a / b)) / c;

N = 0;
while N < x
    N = N + 1;
end

omega_c = 1;

disp(['Analog Butterworth order N = ', num2str(N)]);
disp(['Analog cutoff frequency omega_c = ', num2str(omega_c)]);

w = linspace(0, 2 * omega_s, 500);
Ha = zeros(size(w));

temp = N;
while temp >= 2
    temp = temp - 2;
end

if temp == 0
    N_odd = false;
else
    N_odd = true;
end

if N_odd
    for idx = 1:length(w)
        s = 1j * w(idx);
        s = omega_s / s;
        Hprod = 1;
        for k = 1:(N - 1) / 2
            bk = 2 * sin((2 * k - 1) * pi / (2 * N));
            Hprod = Hprod * (s^2 + bk * omega_c * s + omega_c^2);
        end
        Ha(idx) = abs((omega_c^N) / ((s + omega_c) * Hprod));
    end
else
    for idx = 1:length(w)
        s = 1j * w(idx);
        s = omega_s / s;
        Hprod = 1;
        for k = 1:(N / 2)
            bk = 2 * sin((2 * k - 1) * pi / (2 * N));
            Hprod = Hprod * (s^2 + bk * omega_c * s + omega_c^2);
        end
        Ha(idx) = abs((omega_c^N) / Hprod);
    end
end

Nfft = 1024;
wd = linspace(0, pi, Nfft);
Hd_blt = zeros(size(wd));

for idx = 1:length(wd)
    z = exp(1j * wd(idx));
    s = (2 / T) * (1 - z^-1) / (1 + z^-1);
    s = omega_s / s;
    Hprod = 1;
    for k = 1:N
        thetak = (2 * k + N - 1) * pi / (2 * N);
        ck = omega_c * exp(1j * thetak);
        Hprod = Hprod * (omega_c / (s - ck));
    end
    Hd_blt(idx) = Hprod;
end

figure;
subplot(2, 1, 1);
plot(w, 20 * log10(Ha), 'LineWidth', 1.5);
xlabel('\omega (rad/s)');
ylabel('Magnitude (dB)');
title('Analog Butterworth LPF (Product Form)');
grid on;
ylim([-60 5]);

subplot(2, 1, 2);
fHz = wd * fs / (2 * pi);
plot(fHz, 20 * log10(abs(Hd_blt)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Digital Butterworth LPF (Bilinear Transform)');
grid on;
ylim([-60 5]);
