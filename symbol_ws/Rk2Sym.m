% 定义符号
syms x y theta delta v a u1 u2 L k dt

% 定义状态变量和状态向量
state = [x; y; theta; delta; v; a];

% 定义状态方程
f = [
    v * cos(theta);
    v * sin(theta);
    v * tan(delta) / (L * (1 + k * v^2));
    u1;
    a;
    u2
];

% Runge-Kutta 2 方法离散化
k1 = f;
k2 = subs(f, [x, y, theta, delta, v, a], [x + k1(1) * dt * 0.5, y + k1(2) * dt  * 0.5, theta + k1(3) * dt * 0.5, delta + k1(4) * dt * 0.5, v + k1(5) * dt * 0.5, a + k1(6) * dt * 0.5]);

% 新的状态变量
state_new = state + k2 * dt;

% 简化结果
state_new = simplify(state_new);

% 计算雅各比矩阵
jacobian_matrix = jacobian(state_new, state);

% 显示雅各比矩阵
disp(jacobian_matrix);
