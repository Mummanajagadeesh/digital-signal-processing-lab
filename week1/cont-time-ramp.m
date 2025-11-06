t = -2:0.01:2; 
ramp = t .* (t >= 0); 
plot(t, ramp, 'LineWidth', 2);
xlabel('t');
ylabel('r(t)');
title('Continuous-Time Unit Ramp Signal');
grid on;
axis([-2 2 -0.5 2.5]);