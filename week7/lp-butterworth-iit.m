clc;
clear;

p1 = input('Passband frequency (rad/s): ');
p2 = input('Stopband frequency (rad/s): ');
r1 = input('Passband ripple (dB): ');
r2 = input('Stopband attenuation (dB): ');
fsamp = input('Sampling frequency (Hz): ');

ts = 1 / fsamp;

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
Ha = zeros(size(wset));

odd = false;
tmp = n;
while tmp > 1
    tmp = tmp - 2;
end
if tmp == 1
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
        Ha(i) = abs((wc^n) / ((s + wc) * P));
    else
        k = 1;
        while k <= n / 2
            b = 2 * sin((2 * k - 1) * pi / (2 * n));
            P = P * (s^2 + b * wc * s + wc^2);
            k = k + 1;
        end
        Ha(i) = abs((wc^n) / P);
    end
    i = i + 1;
end

pts = 1024;
wd = linspace(0, pi, pts);

S = zeros(1, n);
k = 0;
while k < n
    t = pi * (2 * k + n + 1) / (2 * n);
    S(k + 1) = wc * exp(1j * t);
    k = k + 1;
end

R = zeros(1, n);
k = 1;
while k <= n
    sk = S(k);
    other = [S(1:k-1) S(k+1:n)];
    den = 1;
    m = 1;
    while m <= length(other)
        den = den * (sk - other(m));
        m = m + 1;
    end
    R(k) = ts * wc^n / den;
    k = k + 1;
end

Hiit = zeros(size(wd));
i = 1;
while i <= length(wd)
    z = exp(1j * wd(i));
    summ = 0;
    k = 1;
    while k <= n
        summ = summ + R(k) / (1 - exp(S(k) * ts) * z^-1);
        k = k + 1;
    end
    Hiit(i) = abs(summ);
    i = i + 1;
end

figure;
subplot(2, 1, 1);
plot(wset, 20 * log10(Ha), 'LineWidth', 1.5);
xlabel('\omega (rad/s)');
ylabel('Magnitude (dB)');
title('Analog Butterworth LPF');
grid on;
ylim([-60 5]);

subplot(2, 1, 2);
faxis = wd * fsamp / (2 * pi);
plot(faxis, 20 * log10(Hiit), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Digital Butterworth LPF (IIT)');
grid on;
ylim([-60 5]);
