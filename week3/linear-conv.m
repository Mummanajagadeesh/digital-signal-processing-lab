x = input('Enter sequence x: ');
x_start = input('Enter starting index of x: ');
h = input('Enter sequence h: ');
h_start = input('Enter starting index of h: ');

% ===== Built-in convolution =====
y_builtin = conv(x, h);

% ===== Manual convolution =====
N1 = length(x);
N2 = length(h);
N = N1 + N2 - 1;
y_manual = zeros(1, N);

for n = 1:N
    for k = 1:N1
        if (n - k + 1 >= 1) && (n - k + 1 <= N2)
            y_manual(n) = y_manual(n) + x(k) * h(n - k + 1);
        end
    end
end

% ===== Index calculation =====
y_start = x_start + h_start;
y_end = y_start + N - 1;
n_y = y_start:y_end;

% ===== Display results =====
disp('Indices:');
disp(n_y);
disp('Built-in conv() Result:');
disp(y_builtin);
disp('Manual Convolution Result:');
disp(y_manual);

% ===== Plot results =====
figure;
subplot(1, 2, 1);
stem(n_y, y_builtin, 'filled');
title('Built-in conv() Result');
xlabel('n');
ylabel('Amplitude');
grid on;

subplot(1, 2, 2);
stem(n_y, y_manual, 'filled');
title('Manual Convolution Result');
xlabel('n');
ylabel('Amplitude');
grid on;
