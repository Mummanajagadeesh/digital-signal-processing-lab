n = 0:20;
% Unit step function
u = ones(size(n)); % u[n] = 1 for n >= 0
% 1. Periodic signal: x[n] = ...,1,2,3,4,...
x1 = repmat([1 2 3 4], 1, ceil(length(n)/4));
x1 = x1(1:length(n));
% 2. Exponential decay: x[n] = 20*(0.9)^n * u[n]
x2 = 20 * (0.9).^n .* u;
% 3. Exponential growth: x[n] = 0.2*(1.2)^n * u[n]
x3 = 0.2 * (1.2).^n .* u;
% 4. Alternating decay: x[n] = (-0.8)^n * u[n]
x4 = (-0.8).^n .* u;
% 5. Negative exponential decay: x[n] = -4*(0.8)^n * u[n]
x5 = -4 * (0.8).^n .* u;
% Plotting
figure;
subplot(3,2,1);
stem(n, x1, 'filled');
title('x[n] = periodic \{1,2,3,4\}');
xlabel('n'); ylabel('x[n]');
grid on;
subplot(3,2,2);
stem(n, x2, 'filled');
title('x[n] = 20(0.9)^n u[n]');
xlabel('n'); ylabel('x[n]');
grid on;
subplot(3,2,3);
stem(n, x3, 'filled');
title('x[n] = 0.2(1.2)^n u[n]');
xlabel('n'); ylabel('x[n]');
grid on;
subplot(3,2,4);
stem(n, x4, 'filled');
title('x[n] = (-0.8)^n u[n]');
xlabel('n'); ylabel('x[n]');
grid on;
subplot(3,2,5);
stem(n, x5, 'filled');
title('x[n] = -4(0.8)^n u[n]');
xlabel('n'); ylabel('x[n]');
grid on;