 = -2:0.01:2; 
u = (t >= 0); 
plot(t, u, 'LineWidth', 2);
xlabel('t');
ylabel('u(t)');
title('Continuous-Time Unit Step Signal');
grid on;
axis([-2 2 -0.2 1.2]);