Xin = input('Enter the frequency-domain sequence: ');
N = length(Xin);

x_dit = (1 / N) * conj(myDIT(conj(Xin)))';
x_dif = (1 / N) * conj(myDIF(conj(Xin)))';
x_sys = ifft(Xin);

m1 = abs(x_dit); p1 = angle(x_dit);
m2 = abs(x_dif); p2 = angle(x_dif);
m3 = abs(x_sys); p3 = angle(x_sys);

n = 0:N-1;

figure;
subplot(3,2,1);
stem(n, m1, 'filled');
title('DIT–IDFT |x[n]|');

subplot(3,2,2);
stem(n, p1, 'filled');
title('DIT–IDFT ∠x[n]');

subplot(3,2,3);
stem(n, m2, 'filled');
title('DIF–IDFT |x[n]|');

subplot(3,2,4);
stem(n, p2, 'filled');
title('DIF–IDFT ∠x[n]');

subplot(3,2,5);
stem(n, m3, 'filled');
title('IFFT |x[n]|');

subplot(3,2,6);
stem(n, p3, 'filled');
title('IFFT ∠x[n]');

disp('IDFT using DIT:');
disp(x_dit.');

disp('IDFT using DIF:');
disp(x_dif.');

disp('Built-in IFFT:');
disp(x_sys);

function Y = myDIT(seq)
L = length(seq);
if L == 1
    Y = seq;
else
    e = seq(1:2:end);
    o = seq(2:2:end);
    E = myDIT(e);
    O = myDIT(o);
    W = exp(-1j * 2 * pi * (0:L/2 - 1) / L);
    Y = [E + W .* O, E - W .* O];
end
end

function Y = myDIF(seq)
L = length(seq);
if L == 1
    Y = seq;
else
    a = seq(1:L/2) + seq(L/2 + 1:L);
    b = (seq(1:L/2) - seq(L/2 + 1:L)) .* exp(-1j * 2 * pi * (0:L/2 - 1) / L);
    A = myDIF(a);
    B = myDIF(b);
    Y = [A, B];
end
end
