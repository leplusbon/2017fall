% Mechanics and Design
% Homework 6
% 2014-10773 ¿Â¿Œ±‘

% for better legibility, used letter J instead of I for moment of inertia

%% problem 1

problem1();

%% problem 2

problem2();

%% (func) problem 1

% problem 1
function problem1()

    E = 10000000;   % Pa
    A = 5;          % m2
    J = 50;         % m4
    L1 = 100;       % m
    L2 = 50;        % m
    L3 = 50;        % m
    
    % elements
    % node1, node2, [cos, sin], E, A, J, L
    elements = [
        1, 2, [1, 0], E, A, J, L1;
        1, 3, [0, -1], E, A, J, L2;
        2, 4, [0, -1], E, A, J, L3
        ];
    
    % boundary conditions
    % first row : index for d being not 0
    % second row : index for known F
    % third row : known F values
    boundary_cond = [
        1, 2, 3, 4, 5, 6;
        1, 2, 3, 4, 5, 6;
        3000, -4000, 0, 0, -3000, 70000
        ];
    
    % obtain global stiffness matrix K_g
    K_g = K_global(4, elements);
    
    % obtain reduced stiffness matrix K_r
    [result, K_r] = solve_for(K_g, boundary_cond);
    
    % results
    d_prob1 = result(:, 1);
    f_prob1 = result(:, 2);
    % buckling analysis : P_cr ~ 10^6 N so it's safe to assume no buckling

    % display results
    disp("-------------------");
    disp("---  Problem 1  ---");
    disp("-------------------");

    disp("displacement of node 1 (m)");
    disp(d_prob1(1:2)');
    disp("roatation of node 1 (deg)");
    disp(d_prob1(3) * 180 / pi);

    disp("displacement of node 2 (m)");
    disp(d_prob1(4:5)');
    disp("rotation of node 2 (deg)");
    disp(d_prob1(6) * 180 / pi);

    disp("force on node 3 (N)");
    disp(f_prob1(7:8)');
    disp("moment on node 3 (N m)");
    disp(f_prob1(9));

    disp("force on node 4 (N)");
    disp(f_prob1(10:11)');
    disp("moment on node 4 (N m)");
    disp(f_prob1(12));
    
end

%% (func) problem 2

% problem 2
function problem2()
    
    E = 200000000000;       % Pa
    b = 0.1;                % m
    h = 0.2;                % m
    A = b * h;              % m2
    J = b * h * h * h / 12; % m4
    L = 1;                  % m
    P = 100000;             % N
    
    % solution for (c) (n = 2)
    res_2 = p2_n(2, E, A, J, L, P);
    dy_2 = res_2(2, 2);

    % solution for (c) (n = 4)
    res_4 = p2_n(4, E, A, J, L, P);
    dy_4 = res_4(2, 3);

    % solution for (c) (n = 10)
    res_10 = p2_n(10, E, A, J, L, P);
    dy_10 = res_10(2, 6);
    
    % solution for (c) (n = 50)
    res_50 = p2_n(50, E, A, J, L, P);
    dy_50 = res_50(2, 26);
    
    % solution for (c) (n = 100)
    res_100 = p2_n(100, E, A, J, L, P);
    dy_100 = res_100(2, 51);
    
    % display results
    disp("-------------------");
    disp("---  Problem 2  ---");
    disp("-------------------");
    
    disp("(c)");
    disp("vertical deflection (m) when n = 2");
    disp(dy_2);
    
    disp("(d)");
    disp("vertical deflection (m) when n = 4");
    disp(dy_4);
    
    disp("vertical deflection (m) when n = 10");
    disp(dy_10);
    
    disp("vertical deflection (m) when n = 50");
    disp(dy_50);
    
    disp("vertical deflection (m) when n = 100");
    disp(dy_100);
    
    % plot the results
    plot(res_2(1, :), res_2(2, :), res_4(1, :), res_4(2, :), ...
        res_10(1, :), res_10(2, :), res_100(1, :), res_100(2, :));
    
end

% generate elements matrix and boundary conditions for problem 2
function deflection = p2_n(n, E, A, J, L, P)

    % number of nodes
    nodes = n + 1;
    
    % elements matrix for n elements
    elements = [(1:n)', (2:(n + 1))', [ones(n, 1), zeros(n, 1)], ...
        E * ones(n, 1), A * ones(n, 1), J * ones(n, 1), ...
        (L / n) * ones(n, 1)];
    
    % boundary conditions for n + 1 nodes, 3 indexes for a node
    unknown_d = [(4:(3 * n + 1)), (3 * n + 3)];
    known_f = [1:(3 * n + 3); zeros(1, 3 * n + 3)];
    
    % concentrated load at x = L/2
    known_f(2, 3 * n / 2 + 2) = -P;
    
    % unknown fx, fy, m at x = 0, unknown fy at x = L;
    known_f(:, 3 * n + 2) = [];
    known_f(:, 1:3) = [];

    boundary_cond = [unknown_d; known_f];
    
    % global stiffness matrix
    K_g = K_global(nodes, elements);
    
    % solve the problem
    [result, K_r] = solve_for(K_g, boundary_cond);
    
    d = result(:, 1);
    
    % deflection in y direction, with x axis
    deflection = [0:(1 / n):1; d(2:3:(3 * n + 2))'];
    
end

%% FEM-solving functions

% obtain the global stiffness matrix
function K_g = K_global(nodes, elements)

    % global stiffness matrix K_global
    % d1, d2, phi for a node
    K_g = zeros(nodes * 3);
    
    % add elements to K_global
    for i = 1:size(elements)
        
        % information about an element
        node1 = elements(i, 1);
        node2 = elements(i, 2);
        c = elements(i, 3);
        s = elements(i, 4);
        E = elements(i, 5);
        A = elements(i, 6);
        J = elements(i, 7);
        L = elements(i, 8);
        
        % some useful factors to reduce # of calculations
        k1 = E * A / L;
        k2 = 12 * E * J / L / L / L;
        k3 = 6 * E * J / L / L;
        k4 = 4 * E * J / L;
        k5 = 2 * E * J / L;
        
        % element stiffness matrix in local coord
        k_element_local = [
            k1, 0, 0, -k1, 0, 0;
            0, k2, k3, 0, -k2, k3;
            0, k3, k4, 0, -k3, k5;
            -k1, 0, 0, k1, 0, 0;
            0, -k2, -k3, 0, k2, -k3;
            0, k3, k5, 0, -k3, k4
            ];
        
        % transformation matrix
        T = [
            c, s, 0, 0, 0, 0;
            -s, c, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            0, 0, 0, c, s, 0;
            0, 0, 0, -s, c, 0;
            0, 0, 0, 0, 0, 1
            ];
        
        % coord transformation from local to global
        k_element_global = T' * k_element_local * T;
        
        % add to K_global
        index = [(3 * node1 - 2):(3 * node1), (3 * node2 - 2):(3 * node2)];
        K_g(index, index, 1) = K_g(index, index, 1) ...
            + k_element_global;
        
    end
    
end

% solve the problem, returns d, f, and K_reduced, reused from HW5
function [result, K_reduced] = solve_for(K_global, boundary_cond)

    % generate the reduced stiffness matrix
    K_reduced = K_global(boundary_cond(1, :), boundary_cond(2, :), 1);
    
    % known F's
    f_reduced = boundary_cond(3, :)';
    
    % obtain unknown d's from known F's
    d_reduced = K_reduced \ f_reduced;
    
    % known d's are all zero
    d = zeros(size(K_global, 2), 1);
    d(boundary_cond(1, :), 1, 1) = d_reduced;
    f = K_global * d;
    
    result = [d, f];
    
end

