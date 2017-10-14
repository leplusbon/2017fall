% Mechanics and Design
% Homework 5
% 2014-10773 ¿Â¿Œ±‘

% 2. 

% number of nodes
nodes = 4;

% some constants
sqrt2 = sqrt(2);
sqrt1_2 = sqrt(0.5);

% parameters / values
E = 1000000000;     % Pa
L = 1;              % m
A = 0.0001;         % m2
F = 1000;           % N

% global stiffness matrix
K_global = zeros(nodes * 2);

% add elements
K_global = add_element(K_global, 1, 4, E * A * sqrt2 / L, ...
    [-sqrt1_2; -sqrt1_2]);
K_global = add_element(K_global, 1, 2, E * A * sqrt2 / L, ...
    [sqrt1_2; sqrt1_2]);
K_global = add_element(K_global, 1, 3, E * A * sqrt2 / L, ...
    [-sqrt1_2; sqrt1_2]);
K_global = add_element(K_global, 3, 4, E * A / L, [0; -1]);
K_global = add_element(K_global, 2, 3, E * A / L, [-1; 0]);

% boundary conditions
% first row : index for d being not 0
% second row : index for known F
% third row : known F values
boundary_cond = [
    1, 2, 3;
    1, 2, 3;
    0, -F, 0
    ];

[result, K_reduced] = solve_for(K_global, boundary_cond);
d = result(:, 1);
f = result(:, 2);

% display the results

% global stiffness matrix
disp("global stiffness matrix (N/m)");
disp(K_global);

% reduced stiffness matrix
disp("reduced stiffness matrix (N/m)");
disp(K_reduced);

% displacement of node 1
disp("displacement of node 1 (m)");
disp(d(1:2));

% displacement of node 2
disp("displacement of node 2 (m)");
disp(d(3:4));

% stress in truss 14
disp("stress in truss 14 (Pa)");
disp([-sqrt1_2; -sqrt1_2]' * (d(7:8) - d(1:2)) * E / L * sqrt2);

% stress in truss 23
disp("stress in truss 23 (Pa)");
disp([-1; 0]' * (d(5:6) - d(3:4)) * E / L);

% add elements between nodes, returns K_global
function K_global = add_element(K_global, node1, node2, k, direction)
    
    % element stiffness matrix in local coord
    K_element_local = [k, -k; -k, k];
    
    % transformation matrix from local coord to global coord
    T = [
        direction(1), direction(2), 0, 0;
        0, 0, direction(1), direction(2)
        ];
    
    % element stiffness matrix in global coord
    K_element_global = T' * K_element_local * T;
    
    index = [2 * node1 - 1, 2 * node1, 2 * node2 - 1, 2 * node2];
    
    % add the element stiffness matrix to K_global
    K_global(index, index, 1) ...
        = K_global(index, index, 1) + K_element_global;
    
end

% solve the problem
function [result, K_reduced] = solve_for(K_global, boundary_cond)
    K_reduced = K_global(boundary_cond(1, :), boundary_cond(2, :), 1);
    f_reduced = boundary_cond(3, :)';
    d_reduced = K_reduced \ f_reduced;
    d = zeros(size(K_global, 2), 1);
    d(boundary_cond(1, :), 1, 1) = d_reduced;
    f = K_global * d;
    result = [d, f];
end

