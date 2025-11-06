n = 0:50;
T = 10; 
W = 4; 
x = mod(n, T) < W; 
stem(n, x, 'filled');
xlabel('n');
ylabel('x[n]');
title('Pulse Train (T=10, W=4)');
grid on;