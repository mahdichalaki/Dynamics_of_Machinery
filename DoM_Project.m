clc; 
clear all;

%% Parameters

% Link lengths
A = 0.5; 
B = 5; 
C = 3.0413813; 
D = 3; 
F = 2; 
G = 3.5; 
H = 2.804015;

% Masses of the links
m2 = 0.5 * 1;
m3 = 5 * 1;
m4 = C * 1;
m5 = 1;
m6 = G * 1;

% Moments of inertia of the links
I2 = (1/12) * m2 * (0.5^2);
I3 = (1/12) * m3 * (5^2);
I4 = (1/12) * m4 * (C^2);
I6 = (1/12) * m6 * (G^2);

% Time, angular velocity, and theta for link AB
t = 0:0.05:10;
ang_speed = 1;
theta = ang_speed * t + pi;

% Fixed positions of points A, D, and E
P1 = [0; 0];    % Point A
P4 = D * [1; 0]; % Point D
P6 = [1.25; 4];  % Point E

%% Static Analysis

% Calculate position of point B (P2)
P2 = A * [cos(theta); sin(theta)]; 

% Calculate intermediate values for point C (P3)
E = sqrt(A^2 + D^2 - 2 * A * D * cos(theta));
alfa = asin(A * sin(theta) ./ E);
beta = acos((E.^2 + C^2 - B^2) ./ (2 * E * C));

% Calculate position of point C (P3)
P3 = [D - C * cos(alfa + beta); C * sin(alfa + beta)];

% Calculate position of point F5 (P5)
P5 = [P2(1,:) + (P3(1,:) - P2(1,:)) * F / B; 
      P2(2,:) + (P3(2,:) - P2(2,:)) * F / B];

% Extract position components
P2_x = P2(1,:); P2_y = P2(2,:); % B
P3_x = P3(1,:); P3_y = P3(2,:); % C
P5_x = P5(1,:); P5_y = P5(2,:); % F5

% Calculate velocities and accelerations of points B, C, and F5
[P2_v, P2_a] = velocity_acceleration(P2_x, P2_y, t);
[P3_v, P3_a] = velocity_acceleration(P3_x, P3_y, t);
[P5_v, P5_a] = velocity_acceleration(P5_x, P5_y, t);

% Angle between EF and vertical
gamma = atan((P5(1,:) - P6(1)) ./ (P6(2,:) - P5(2,:)));
gamma_dot = diff(gamma) ./ diff(t);
gamma_double_dot = diff(gamma_dot) ./ diff(t(1:end-1));

% Distance EF
EF = sqrt((P5(1,1:end-2) - P6(1)).^2 + (P5(2,1:end-2) - P6(2)).^2);

% Velocity and acceleration of point F6
P8_v = abs(sqrt((P5(1,end-1) - P6(1)).^2 + (P5(2,end-1) - P6(2)).^2) .* gamma_dot);
P8_a = sqrt((EF .* (gamma_dot(1:end-1).^2)).^2 + (EF .* gamma_double_dot).^2);

%% Animation
for i = 1:length(t)
    subplot(2,1,1);
    hold on;
    
    % Draw circles representing link joints
    viscircles(P1', 0.05); 
    viscircles(P2(:,i)', 0.05); 
    viscircles(P3(:,i)', 0.05); 
    viscircles(P4', 0.05); 
    viscircles(P5(:,i)', 0.05); 
    viscircles(P6', 0.05);

    % Draw connecting bars
    line([P1(1) P2(1,i)], [P1(2) P2(2,i)], 'Color', 'b'); % Link AB
    line([P2(1,i) P3(1,i)], [P2(2,i) P3(2,i)], 'Color', 'r'); % Link BC
    line([P3(1,i) P4(1)], [P3(2,i) P4(2)], 'Color', 'g'); % Link CD
    line([P6(1) P5(1,i)], [P6(2) P5(2,i)], 'Color', 'k'); % Link EF

    axis equal;
    xlim([-5 8]); ylim([-2 7]);

    % Display time
    text(-2, 6, ['Time elapsed: ', num2str(t(i)), ' s']);
    
    pause(0.005);
    clf;
end

%% Plot velocities and accelerations
plot_velocity_acceleration(t, P2_v, P3_v, P5_v, P8_v, 'Speed');
plot_velocity_acceleration(t, P2_a, P3_a, P5_a, P8_a, 'Acceleration');

%% Dynamics Analysis

% Forces and Torque Calculation
[forces, torque] = calculate_dynamics(t, P2_a, P3_a, P5_a, gamma, gamma_dot, gamma_double_dot, m2, m3, m4, m5, m6, I2, I3, I4, I6);

% Plot Forces and Torque
plot_forces_torque(t, forces, torque);

%% Function Definitions

function [velocity, acceleration] = velocity_acceleration(x, y, t)
    vx = diff(x) ./ diff(t);
    vy = diff(y) ./ diff(t);
    velocity = sqrt(vx.^2 + vy.^2);
    ax = diff(vx) ./ diff(t(1:end-1));
    ay = diff(vy) ./ diff(t(1:end-1));
    acceleration = sqrt(ax.^2 + ay.^2);
end

function plot_velocity_acceleration(t, v1, v2, v3, v4, title_str)
    figure;
    subplot(2,2,1); plot(t(1:end-1), v1); title([title_str, ' of B']);
    subplot(2,2,2); plot(t(1:end-1), v2); title([title_str, ' of C']);
    subplot(2,2,3); plot(t(1:end-1), v3); title([title_str, ' of F5']);
    subplot(2,2,4); plot(t(1:end-1), v4); title([title_str, ' of F6']);
    grid on;
end

function [forces, torque] = calculate_dynamics(t, a2, a3, a5, gamma, gamma_dot, gamma_double_dot, m2, m3, m4, m5, m6, I2, I3, I4, I6)
    alpha3 = diff(gamma) ./ diff(t);
    alpha4 = diff(alpha3) ./ diff(t(1:end-1));
    forces = m2 .* a2 + m3 .* a3 + m5 .* a5;
    torque = I2 * 0 + I3 * alpha3 + I4 * alpha4 + I6 * gamma_double_dot;
end

function plot_forces_torque(t, forces, torque)
    figure;
    subplot(1,2,1); plot(t(1:end-2), forces); title('Forces');
    subplot(1,2,2); plot(t(1:end-2), torque); title('Torque');
    grid on;
end

