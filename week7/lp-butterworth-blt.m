clc;
clear;

p1 = input('Passband frequency (rad/s): ');
p2 = input('Stopband frequency (rad/s): ');
r1 = input('Passband ripple (dB): ');
r2 = input('Stopband attenuation (dB): ');
fsamp = input('Sampling frequency (Hz): ');

ts = 1 / fsamp;

p1 = (2 / ts) * tan((p1 * ts) / 2);
p2 = (2 / ts) * tan((p2 * ts) / 2);

dp = 10^(-r1 / 20);
ds = 10^(-r2 / 20);

A = (1 / ds^2) - 1;
B = (1 / dp^2) - 1;
C = log10(p2 / p1);
val = log10(sqrt(A / B)) / C;

n = 1;
while n < val
    n = n + 1;
end

wc = p1 / ((1 / dp^2 - 1)^(1 / (2 * n)));

disp(['N = ', num2str(n)]);
disp(['wc = ', num2str(wc)]);

wset = linspace(0, 2 * p2, 500);
Hval = zeros(size(wset));

odd = false;
temp = n;
while temp > 1
    temp = temp - 2;
end
if temp == 1
    odd = true;
end

i = 1;
while i <= length(wset)
    s = 1j * wset(i);
    P = 1;
    if odd
        k = 1;
        while k <= (n - 1) / 2
            b = 2 * sin((2 * k - 1) * pi / (2 * n));
            P = P * (s^2 + b * wc * s + wc^2);
            k = k + 1;
        end
        Hval(i) = abs((wc^n) / ((s + wc) * P));
    else
        k = 1;
        while k <= n / 2
            b = 2 * sin((2 * k - 1) * pi / (2 * n));
            P = P * (s^2 + b * wc * s + wc^2);
            k = k + 1;
        end
        Hval(i) = abs((wc^n) / P);
    end
    i = i + 1;
end

pts = 1024;
wd = linspace(0, pi, pts);
Hdb = zeros(size(wd));

j = 1;
while j <= length(wd)
    z = exp(1j * wd(j));
    s = (2 / ts) * ((1 - z^-1) / (1 + z^-1));
    P = 1;
    k = 1;
    while k <= n
        t = ((2 * k) + n - 1) * pi / (2 * n);
        ck = wc * exp(1j * t);
        P = P * (wc / (s - ck));
        k = k + 1;
    end
    Hdb(j) = P;
    j = j + 1;
end

figure;
subplot(2, 1, 1);
plot(wset, 20 * log10(Hval), 'LineWidth', 1.5);
xlabel('\omega (rad/s)');
ylabel('Magnitude (dB)');
title('Analog Butterworth LPF');
grid on;
ylim([-60 5]);

subplot(2, 1, 2);
faxis = wd * fsamp / (2 * pi);
plot(faxis, 20 * log10(abs(Hdb)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Digital Butterworth LPF (BLT)');
grid on;
ylim([-60 5]);
