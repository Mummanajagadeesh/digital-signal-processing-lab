x = input('Enter input sequence x (row vector): ');
n0x = input('Enter the index of x[1] (first element of x): ');
h = input('Enter impulse response h (row vector): ');
n0h = input('Enter the index of h[1] (first element of h): ');
L = input('Enter block length L: ');

Nx = length(x);
Nh = length(h);
M = L + Nh - 1;

numblocks = 0;
while numblocks * L < Nx
    numblocks = numblocks + 1;
end

xpadded = [x, zeros(1, numblocks * L - Nx)];
y = zeros(1, Nx + Nh - 1);

for b = 0:numblocks - 1
    startidx = b * L + 1;
    endidx = startidx + L - 1;
    xblock = xpadded(startidx:endidx);
    yblock = zeros(1, M);
    for n = 1:M
        s = 0;
        for k = 1:Nh
            if (n - k + 1 > 0 && n - k + 1 <= length(xblock))
                s = s + h(k) * xblock(n - k + 1);
            end
        end
        yblock(n) = s;
    end
    ystart = startidx;
    yend = min(startidx + M - 1, length(y));
    blkend = yend - startidx + 1;
    y(ystart:yend) = y(ystart:yend) + yblock(1:blkend);
end

n0y = n0x + n0h;
ny = n0y:n0y + length(y) - 1;

yd = zeros(1, Nx + Nh - 1);
for n = 1:(Nx + Nh - 1)
    s = 0;
    for k = 1:Nh
        m = n - k + 1;
        if (m > 0 && m <= Nx)
            s = s + h(k) * x(m);
        end
    end
    yd(n) = s;
end

nyd = n0x + n0h : n0x + n0h + length(yd) - 1;

disp('Result of convolution using Overlap-Add:');
disp(['Indices: ', mat2str(ny)]);
disp(['Values: ', mat2str(y)]);
disp('Result of convolution using Direct Method:');
disp(['Indices: ', mat2str(nyd)]);
disp(['Values: ', mat2str(yd)]);

subplot(4,1,1);
stem(n0x:n0x + Nx - 1, x, 'filled', 'LineWidth', 2);
title('Input Sequence x[n]');
xlabel('n'); ylabel('x[n]');
set(gca, 'FontWeight', 'bold');

subplot(4,1,2);
stem(n0h:n0h + Nh - 1, h, 'filled', 'LineWidth', 2);
title('Impulse Response h[n]');
xlabel('n'); ylabel('h[n]');
set(gca, 'FontWeight', 'bold');

subplot(4,1,3);
stem(ny, y, 'filled', 'LineWidth', 2);
title('Output y[n] using Overlap-Add');
xlabel('n'); ylabel('y[n]');
set(gca, 'FontWeight', 'bold');

subplot(4,1,4);
stem(nyd, yd, 'filled', 'LineWidth', 2);
title('Output y[n] using Direct Convolution');
xlabel('n'); ylabel('y[n]');
set(gca, 'FontWeight', 'bold');
