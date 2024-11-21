%Code by Jerry Paul Varghese, Adcom. Mat.N. 453553
%Dated: 21/11/2024
% Parameters
m2 = 1;  % Mass at node 2
m3 = 1;  % Mass at node 3
k1 = 10^7;  % Stiffness between node 1 and node 2
k2 = 1;    % Stiffness between node 2 and node 3
omega_p = 1.2;  % Prescribed frequency for displacement at node 1
dt = 0.1;  % Time step size
total_time = 10;  %Total simulation time
time = 0:dt:total_time;  % Time vector

% Newmark parameters for trapezoidal method
beta = 0.25;
gamma = 0.5;

% Mass and stiffness matrices
M = [m2, 0;
     0, m3];
 
K = [k1 + k2, -k2;
    -k2, k2];

% Initial conditions
u = zeros(2, length(time));  % Displacement at nodes 2 and 3
udot = zeros(2, length(time));  % Velocity at nodes 2 and 3
uddot = zeros(2, length(time));  % Acceleration at nodes 2 and 3

% Initial acceleration
F_ext = [0; 0];  % No external forces at nodes 2 and 3
uddot(:, 1) = M \ (F_ext - K * u(:, 1));

% Time stepping using trapezoidal method
for i = 2:length(time)
    % Prescribed displacement at node 1 (inverted sine wave)
    u1 = -sin(omega_p * time(i));
    
    % Effective force vector
    F_eff = F_ext - [k1 * u1; 0] + M * ((1 / (beta * dt^2)) * u(:, i-1) + (1 / (beta * dt)) * udot(:, i-1) + (1 / (2 * beta) - 1) * uddot(:, i-1));
    
    % Effective stiffness matrix
    K_eff = K + (1 / (beta * dt^2)) * M;
    
    % Solve for displacement at nodes 2 and 3
    u(:, i) = K_eff \ F_eff;
    
    % Calculate velocity and acceleration
    udot(:, i) = (gamma / (beta * dt)) * (u(:, i) - u(:, i-1)) - (gamma / beta - 1) * udot(:, i-1) - (gamma / (2 * beta) - 1) * dt * uddot(:, i-1);
    uddot(1, i) = (1 / (beta * dt^2)) * (u(1, i) - u(1, i-1)) - (1 / (beta * dt)) * udot(1, i-1) - (1 / (2 * beta) - 1) * uddot(1, i-1);  % Node 2 using Newmark's force balance
    uddot(2, i) = (udot(2, i) - udot(2, i-1)) / dt;  % Node 3 using backward difference of velocity
end

% Plotting results in a single figure
figure;
sgtitle('Newmarks Method');


% Displacement at nodes 2 and 3
subplot(2, 2, 1);
plot(time, u(1, :), 'r', 'DisplayName', 'Displacement at Node 2');
hold on;
plot(time, u(2, :), 'b', 'DisplayName', 'Displacement at Node 3');
xlabel('Time (s)');
ylabel('Displacement (m)');
legend;
grid on;
title('Displacement');


% Velocity at nodes 2 and 3
subplot(2, 2, 2);
plot(time, udot(1, :), 'r', 'DisplayName', 'Velocity at Node 2');
hold on;
plot(time, udot(2, :), 'b', 'DisplayName', 'Velocity at Node 3');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend;
grid on;
title('Velocity');

% Acceleration at Node 2
subplot(2, 2, 3);
plot(time, uddot(1, :), 'r', 'DisplayName', 'Acceleration at Node 2');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend;
grid on;
title('Acceleration');

% Acceleration at Node 3
subplot(2, 2, 4);
plot(time, uddot(2, :), 'r', 'DisplayName', 'Acceleration at Node 3 ');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend;
grid on;
title('Acceleration');
xlim([0, total_time]);
ylim([-5, 5]);  % Adjusted limits for acceleration similar to paper's reference plot