n = -10:10; 
u = (n >= 0); 
stem(n, u, 'filled');
xlabel('n');
ylabel('u[n]');
title('Discrete Unit Step Signal');
grid on;
