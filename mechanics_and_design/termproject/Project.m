

%% given values

% thickness (m)
thickness = 100;

% Young's modulus (Pa)
E = 200e+9;

% Poisson's ratio
nu = 0.3;

% weight per unit length (N/m)
w = 100e+3;

% stress = D * strain, plane strain
D = [
    1-nu, nu, 0;
    nu, 1-nu, 0;
    0, 0, (1-2*nu)/2
    ] * E / (1 + nu) / (1 - 2 * nu);

% boundary (m)
boundary = [
    0, 0;
    20, 0;
    20, 1;
    0, 1
    ];

% given nodes
initial_node = zeros(42, 2);
% for i = 1:21
%     initial_node(2*i-1, :) = [i-1, 0];
%     initial_node(2*i, :) = [i-1, 1];
% end

for i = 1:21
    initial_node(i, :) = [i-1, 0];
    initial_node(i+21, :) = [i-1, 1];
end

% given elements from given nodes
initial_element = zeros(40, 3);
% for i = 1:20
%     initial_element(2*i-1, :) = [2*i-1, 2*i+1, 2*i];
%     initial_element(2*i, :) = [2*i, 2*i+1, 2*i+2];
% end

initial_element = delaunay(initial_node);

triplot(initial_element, initial_node(:, 1), initial_node(:, 2), 'k');
labels = num2str((1:42)', '%d');
text(initial_node(:, 1), initial_node(:, 2), labels, 'horizontal','left', 'vertical','bottom');
%axis [-5 6 0 11]
axis equal

K_global = global_stiffness(initial_node, initial_element, thickness, D);

K_reduced = K_global(5:80, 5:80);

calc(K_reduced);

function calc(matrix)
    tic
    for i = 1:100000
        k = rand(76, 1);
        f = matrix \ k;
    end
    t = toc;
    disp(t);
end


% obtain element stiffness matrix for triangular element
function K_element = element_stiffness(nodes, thickness, D)

    area = polyarea(nodes(:, 1), nodes(:, 2));
    b1 = nodes(2, 2) - nodes(3, 2);
    b2 = nodes(3, 2) - nodes(1, 2);
    b3 = nodes(1, 2) - nodes(2, 2);
    c1 = nodes(3, 1) - nodes(2, 1);
    c2 = nodes(1, 1) - nodes(3, 1);
    c3 = nodes(2, 1) - nodes(1, 1);
    B = [
        b1, 0 , b2, 0 , b3, 0 ;
        0 , c1, 0 , c2, 0 , c3;
        c1, b1, c2, b2, c3, b3
    ] / (2 * area);
    K_element = thickness * area * B' * D * B;
    
end

% obtain global stiffness matrix for triangular mesh
function K_global = global_stiffness(nodes, elements, thickness, D)

    K_global = zeros(2*size(nodes, 1));

    %iteration for each element
    for element = elements'
        
        element = element';
        
        % global stiffness matrix index corresponding to the element
        index = [2*element-1; 2*element];
        index = index(:)';
        
        % element stiffness matrix for a element
        K_element = element_stiffness(nodes(element, :), thickness, D);
        
        % add element stiffness matrix to the global stiffness matrix
        K_global(index, index) = K_global(index, index) + ...
            K_element;
        
    end
    
end

