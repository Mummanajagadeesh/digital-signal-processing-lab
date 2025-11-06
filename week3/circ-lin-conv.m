x = input('Enter sequence x: ');
x_start = input('Enter starting index of x: ');
h = input('Enter sequence h: ');
h_start = input('Enter starting index of h: ');

N1 = length(x);
N2 = length(h);
N_lin = N1 + N2 - 1;

x_pad = [x, zeros(1, N_lin - N1)];
h_pad = [h, zeros(1, N_lin - N2)];

% ===== Linear convolution via circular convolution =====
y_circ_for_linear = zeros(1, N_lin);
for n = 1:N_lin
    for m = 1:N_lin
        k = n - m;
        while k < 0
            k = k + N_lin;
        end
        while k >= N_lin
            k = k - N_lin;
        end
        y_circ_for_linear(n) = y_circ_for_linear(n) + x_pad(m) * h_pad(k + 1);
    end
end

% ===== Manual linear convolution =====
y_linear = zeros(1, N_lin);
for n = 1:N_lin
    for k = 1:N1
        if (n - k + 1 >= 1) && (n - k + 1 <= N2)
            y_linear(n) = y_linear(n) + x(k) * h(n - k + 1);
        end
    end
end

n_axis = (x_start + h_start) : (x_start + h_start + N_lin - 1);

disp('Indices:');
disp(n_axis);
disp('Manual Linear Convolution:');
disp(y_linear);
disp('Linear Convolution via Circular Convolution + Padding:');
disp(y_circ_for_linear);

figure;
subplot(1, 2, 1);
stem(n_axis, y_linear, 'filled', 'LineWidth', 1.5);
title('Manual Linear Convolution');
xlabel('n'); ylabel('Amplitude');
grid on;

subplot(1, 2, 2);
stem(n_axis, y_circ_for_linear, 'filled', 'LineWidth', 1.5);
title('Linear Convolution via Circular Convolution + Padding');
xlabel('n'); ylabel('Amplitude');
grid on;
