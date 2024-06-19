clear, clc;
syms x y phi u v omega T_s a delta k_f k_r l_f l_r m I_z

% Define state and input vectors
X = [x; y; phi; u; v; omega];
U = [a; delta];

% Define the function F
F = [
    x + T_s * (u * cos(phi) - v * sin(phi));
    y + T_s * (v * cos(phi) + u * sin(phi));
    phi + T_s * omega;
    u + T_s * a;
    (m * u * v + T_s * (l_f * k_f - l_r * k_r) * omega - T_s * k_f * delta * u - T_s * m * u^2 * omega) / (m * u - T_s * (k_f + k_r));
    (I_z * u * omega + T_s * (l_f * k_f - l_r * k_r) * v - T_s * l_f * k_f * delta * u) / (I_z * u - T_s * (l_f^2 * k_f + l_r^2 * k_r))
];

% Compute Jacobian matrices
J_X = jacobian(F, X);  % Jacobian with respect to state variables
J_U = jacobian(F, U);  % Jacobian with respect to input variables

% Display the Jacobian matrices
disp('Jacobian with respect to state variables (X):');
disp(J_X);
disp('Jacobian with respect to input variables (U):');
disp(J_U);
