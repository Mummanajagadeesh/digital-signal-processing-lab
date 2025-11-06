n = -20:20;
N = 4; 
x = mod(n, N) == 0; 
stem(n, x, 'filled');
xlabel('n');
ylabel('x[n]');
title('Impulse Train with Period N = 4');
grid on;