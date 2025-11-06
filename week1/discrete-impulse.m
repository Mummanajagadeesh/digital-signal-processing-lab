n = -10:10;
delta = (n == 0);
stem(n, delta, 'filled');
xlabel('n');
ylabel('\delta[n]');
title('Discrete Impulse Signal');
grid on;