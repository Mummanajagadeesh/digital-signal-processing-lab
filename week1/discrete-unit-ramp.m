n = -10:10;
ramp = n .* (n >= 0); 
stem(n, ramp, 'filled');
xlabel('n');
ylabel('r[n]');
title('Discrete Unit Ramp Signal');
grid on;
