clc;
clear;

omega_p = input('Enter passband edge frequency (rad/s): ');
omega_s = input('Enter stopband edge frequency (rad/s): ');
Ap = input('Enter passband ripple (dB): ');
As = input('Enter stopband attenuation (dB): ');
fs = input('Enter sampling frequency (Hz): ');

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
eps = sqrt(b);
c = acosh(omega_s / omega_p);
x = acosh(sqrt(a / b)) / c;

N = 0;
while N < x
    N = N + 1;
end

disp(['Analog Chebyshev Filter order N = ', num2str(N)]);

temp = N;
while temp >= 2
    temp = temp - 2;
end

if temp == 0
    N_odd = false;
else
    N_odd = true;
end

omega_p = 1;
yN = ((sqrt(1 + 1 / eps^2) + 1 / eps)^(1 / N) - (sqrt(1 + 1 / eps^2) + 1 / eps)^(-1 / N)) / 2;

w = linspace(0, 2 * omega_s, 500);
Ha = zeros(size(w));

if N_odd
    for idx = 1:length(w)
        s = 1j * w(idx);
        s = omega_s / s;
        H = (omega_p^N * yN) / (s + omega_p * yN);
        num = 1;
        deno = 1;
        for k = 1:(N - 1) / 2
            ck = yN^2 + (cos((2 * k - 1) * pi / (2 * N)))^2;
            bk = 2 * yN * sin((2 * k - 1) * pi / (2 * N));
            num = num * ck;
            deno = deno * (s^2 + bk * omega_p * s + ck * omega_p^2);
        end
        Ha(idx) = abs(H * num / deno);
    end
else
    for idx = 1:length(w)
        s = 1j * w(idx);
        s = omega_s / s;
        H = omega_p^N / sqrt(1 + eps^2);
        num = 1;
        deno = 1;
        for k = 1:(N / 2)
            ck = yN^2 + (cos((2 * k - 1) * pi / (2 * N)))^2;
            bk = 2 * yN * sin((2 * k - 1) * pi / (2 * N));
            num = num * ck;
            deno = deno * (s^2 + bk * omega_p * s + ck * omega_p^2);
        end
        Ha(idx) = abs(H * num / deno);
    end
end

Nfft = 1024;
wd = linspace(0, pi, Nfft);
Hd_blt = zeros(size(wd));

for idx = 1:length(wd)
    z = exp(1j * wd(idx));
    s = (2 / T) * (1 - z^-1) / (1 + z^-1);
    s = omega_s / s;
    if N_odd
        H = (omega_p^N * yN) / (s + omega_p * yN);
        num = 1;
        deno = 1;
        for k = 1:(N - 1) / 2
            ck = yN^2 + (cos((2 * k - 1) * pi / (2 * N)))^2;
            bk = 2 * yN * sin((2 * k - 1) * pi / (2 * N));
            num = num * ck;
            deno = deno * (s^2 + bk * omega_p * s + ck * omega_p^2);
        end
        Hd_blt(idx) = H * num / deno;
    else
        H = omega_p^N / sqrt(1 + eps^2);
        num = 1;
        deno = 1;
        for k = 1:(N / 2)
            ck = yN^2 + (cos((2 * k - 1) * pi / (2 * N)))^2;
            bk = 2 * yN * sin((2 * k - 1) * pi / (2 * N));
            num = num * ck;
            deno = deno * (s^2 + bk * omega_p * s + ck * omega_p^2);
        end
        Hd_blt(idx) = H * num / deno;
    end
end

figure;
subplot(2, 1, 1);
plot(w, 20 * log10(Ha), 'LineWidth', 1.5);
xlabel('\omega (rad/s)');
ylabel('Magnitude (dB)');
title('Analog Chebyshev Type-I HPF');
grid on;
ylim([-60 5]);

subplot(2, 1, 2);
fHz = wd * fs / (2 * pi);
plot(fHz, 20 * log10(abs(Hd_blt)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Digital Chebyshev Type-I HPF (Bilinear Transform)');
grid on;
ylim([-60 5]);
