% Mechanics and Design
% Homework 6
% 2014-10773 ¿Â¿Œ±‘

function problem1()
    E = 10000000;   % Pa
    A = 5;          % m2
    J = 50;         % m4
    L1 = 100;       % m
    L2 = 50;        % m
    L3 = 50;        % m
    elements = [
        [1, 2, [1; 0], E, A, J, L1];
        [1, 3, [0; -1], E, A, J, L2];
        [2, 4, [0; -1], E, A, J, L3]
        ];
    K = K_global(4, elements);
end

function K_global = K_global(nodes, elements)
    % global stiffness matrix K_global
    K_global = zeros(nodes * 3);
end

