% Mechanics and Design
% Homework 5
% 2014-10773 ¿Â¿Œ±‘

% 1. 

% number of nodes
nodes = 4;

% some constants
k = 100;
sqrt1_2 = sqrt(0.5);

% global stiffness matrix
K_global = zeros(nodes * 2);

% add elements
K_global = add_element(K_global, 1, 4, k, [sqrt1_2; -sqrt1_2]);     % AD
K_global = add_element(K_global, 2, 4, k, [0, -1]);                 % BD
K_global = add_element(K_global, 3, 4, k, [-sqrt1_2; -sqrt1_2]);    % CD



% add elements between nodes, returns K_global
function K_global = add_element(K_global, node1, node2, k, direction)
    
    % noralize direction vector
    %norm2 = direction' * direction;
    %direction = direction / sqrt(norm2);
    
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

% generate reduced stiffness matrix
function K_reduced = solve_for(boundary_cond)
    boundary_cond(1)
end