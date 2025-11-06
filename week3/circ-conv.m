x = input('Enter sequence x: ');
x_start = input('Enter starting index of x: ');
h = input('Enter sequence h: ');
h_start = input('Enter starting index of h: ');

N1 = length(x);
N2 = length(h);
N = max(N1, N2);

if N1 < N
    x = [x, zeros(1, N - N1)];
end
if N2 < N
    h = [h, zeros(1, N - N2)];
end

% Built-in circular convolution
y_builtin = cconv(x, h, N);

% Manual circular convolution
y_manual = zeros(1, N);
for n = 1:N
    for k = 1:N
        t = n - k;
        while t < 0
            t = t + N;
        end
        while t >= N
            t = t - N;
        end
        y_manual(n) = y_manual(n) + x(k) * h(t + 1);
    end
end

n_y = 0:N-1;

disp('Indices:');
disp(n_y);
disp('Built-in Circular Convolution:');
disp(y_builtin);
disp('Manual Circular Convolution:');
disp(y_manual);

figure;
subplot(1, 2, 1);
stem(n_y, y_builtin, 'filled');
title('Built-in Circular Convolution (cconv)');
xlabel('n');
ylabel('Amplitude');
grid on;

subplot(1, 2, 2);
stem(n_y, y_manual, 'filled');
title('Manual Circular Convolution');
xlabel('n');
ylabel('Amplitude');
grid on;

% ===== Linear Convolution using Circular Convolution with Padding =====
N_lin = N1 + N2 - 1;
x_padded = [x, zeros(1, N_lin - N1)];
h_padded = [h, zeros(1, N_lin - N2)];

y_circ_for_linear = cconv(x_padded, h_padded, N_lin);
y_linear = conv(x, h);
n_axis = 0:N_lin - 1;

disp('Linear Convolution via Circular Convolution + Padding:');
disp(y_circ_for_linear);

figure;
subplot(1, 2, 1);
stem(n_axis, y_linear, 'filled', 'LineWidth', 1.5);
title('Manual Linear Convolution');
xlabel('n');
ylabel('Amplitude');
grid on;

subplot(1, 2, 2);
stem(n_axis, y_circ_for_linear, 'filled', 'LineWidth', 1.5);
title('Linear Convolution via Circular Convolution + Padding');
xlabel('n');
ylabel('Amplitude');
grid on;
