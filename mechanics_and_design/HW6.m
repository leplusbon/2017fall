% Mechanics and Design
% Homework 6
% 2014-10773 ¿Â¿Œ±‘

%% problem 1

result_prob1 = problem1();
d_prob1 = result_prob1(:, 1);
f_prob1 = result_prob1(:, 2);
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

%% problem 2

%% functions

% problem 1
function result = problem1()
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
    
    boundary_cond = [
        1, 2, 3, 4, 5, 6;
        1, 2, 3, 4, 5, 6;
        3000, -4000, 0, 0, -3000, 70000
        ];
    
    % obtain global stiffness matrix K_g
    K_g = K_global(4, elements);
    
    % obtain reduced stiffness matrix K_r
    [result, K_r] = solve_for(K_g, boundary_cond);
end

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
