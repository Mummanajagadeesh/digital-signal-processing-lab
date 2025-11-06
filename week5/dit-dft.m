y = input('Enter the input sequence: ');
n = length(y);

Ydit = myDIT(y);
Ydif = myDIF(y);
Ysys = fft(y);

disp('DIT FFT Result:'); disp(Ydit.');
disp('DIF FFT Result:'); disp(Ydif.');
disp('Inbuilt FFT Result:'); disp(Ysys.');

k = 0:n-1;

figure;
subplot(3,2,1);
stem(k, abs(Ydit), 'filled');
title('DIT |X[k]|');

subplot(3,2,2);
stem(k, angle(Ydit), 'filled');
title('DIT ∠X[k]');

subplot(3,2,3);
stem(k, abs(Ydif), 'filled');
title('DIF |X[k]|');

subplot(3,2,4);
stem(k, angle(Ydif), 'filled');
title('DIF ∠X[k]');

subplot(3,2,5);
stem(k, abs(Ysys), 'filled');
title('Built-in |X[k]|');

subplot(3,2,6);
stem(k, angle(Ysys), 'filled');
title('Built-in ∠X[k]');

function Y = myDIT(seq)
m = length(seq);
if m == 1
    Y = seq;
else
    e = seq(1:2:end);
    o = seq(2:2:end);
    E = myDIT(e);
    O = myDIT(o);
    W = exp(-1j * 2 * pi * (0:m/2 - 1) / m);
    Y = [E + W .* O, E - W .* O];
end
end

function Y = myDIF(seq)
seq = seq(:).';
m = length(seq);
if floor(log2(m)) ~= log2(m)
    error('Size must be power of 2');
end
Ytemp = difRec(seq);
idx = 0:m-1;
r = revBits(idx, m);
Y = zeros(1, m);
for p = 1:m
    Y(r(p) + 1) = Ytemp(p);
end
end

function Y = difRec(s)
m = length(s);
if m == 1
    Y = s;
else
    a = s(1:m/2) + s(m/2+1:m);
    b = (s(1:m/2) - s(m/2+1:m)) .* exp(-1j * 2 * pi * (0:m/2 - 1) / m);
    A = difRec(a);
    B = difRec(b);
    Y = [A, B];
end
end

function r = revBits(idx, m)
bits = log2(m);
r = zeros(size(idx));
for q = 1:length(idx)
    val = idx(q);
    rev = 0;
    for j = 1:bits
        b = bitand(val, 1);
        rev = bitor(bitshift(rev, 1), b);
        val = bitshift(val, -1);
    end
    r(q) = rev;
end
end
