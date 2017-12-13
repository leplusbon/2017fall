

% given values

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
for i = 1:21
    initial_node(2*i-1, :) = [i-1, 0];
    initial_node(2*i, :) = [i-1, 1];
end

% given elements from given nodes
initial_element = zeros(40, 3);
for i = 1:20
    initial_element(2*i-1, :) = [2*i-1, 2*i+1, 2*i];
    initial_element(2*i, :) = [2*i, 2*i+1, 2*i+2];
end

figure(1);
triplot(initial_element, initial_node(:, 1), initial_node(:, 2), 'k');
labels = num2str((1:42)', '%d');
text(initial_node(:, 1), initial_node(:, 2), labels, 'horizontal','left', 'vertical','bottom');
%axis [-5 6 0 11]
axis equal

% global stiffness matrix
K_global = global_stiffness(initial_node, initial_element, thickness, D);

% boundary conditions given
% first row : indexes for unknown displacement
% second row : forces corresponding to the index
initial_boundary_cond = zeros(2, 76);
initial_boundary_cond(1, :) = 5:80;
for i = 4:4:size(initial_boundary_cond, 2)
    initial_boundary_cond(2, i) = -100000;
end

% solve for given mesh
[d, f] = solve_for(K_global, initial_boundary_cond);

d = reshape(d, 2, [])';
dx = d(:, 1);
dy = d(:, 2);

f = reshape(f, 2, [])';
fx = f(:, 1);
fy = f(:, 2);

% plot deformed shape with alpha=500000
deformed_shape = initial_node + 500000 * d;
figure(2);
triplot(initial_element, deformed_shape(:, 1), deformed_shape(:, 2), 'k');
xlabel("x (m)");
ylabel("y (m)");
axis equal

% finding the node number with maximum displacement
displacement_magnitude = sqrt(dx .* dx + dy .* dy);
[maxd, maxd_node] = max(displacement_magnitude);
disp(dx(maxd_node));
disp(dy(maxd_node));

% stress of each element
st = element_stress(initial_node, initial_element, dx, dy, D);

% stress of each node
st = node_stress(initial_node, initial_element, st);

% von mises stress
s_v = von_mises(st, nu);
figure(1);
trisurf(initial_element, initial_node(:, 1), initial_node(:, 2), s_v);

% tresca stress
s_t = tresca(st, nu);
figure(2);
trisurf(initial_element, initial_node(:, 1), initial_node(:, 2), s_t);

% maximum values
max_von_mises = max(s_v);
max_tresca = max(s_t);

res=boundary_cond(2, w);
elements=res{1,1};
nodes=res{1, 2};

node_density = [30];

for i = 1:size(node_density, 2)
    res = boundary_cond(node_density(1, i), w);
    elements = res{1, 1};
    nodes = res{1, 2};
    bcond = res{1, 3};
    K_global = global_stiffness(nodes, elements, thickness, D);
    disp("ready to solve for node density(#/m) : ");
    disp(node_density(1, i));
    [d, f] = solve_for(K_global, bcond);
    disp("solved for node density(#/m) : ");
    disp(node_density(1, i));
    d = reshape(d, 2, [])';
    dx = d(:, 1);
    dy = d(:, 2);

    f = reshape(f, 2, [])';
    fx = f(:, 1);
    fy = f(:, 2);

    stress = element_stress(nodes, elements, dx, dy, D);
    stress = node_stress(nodes, elements, stress);
    s_v = von_mises(stress, nu);
    s_t = tresca(stress, nu);
    [max_d(i), node_num_d(i)] = max(sqrt(dx.^2 + dy.^2));
    [max_s_v(i), node_num_s_v(i)] = max(s_v);
    [max_s_t(i), node_num_s_t(i)] = max(s_t);
    loc_d(i, :) = nodes(node_num_d(i), :);
    loc_s_v(i, :) = nodes(node_num_s_v(i), :);
    loc_s_t(i, :) = nodes(node_num_s_t(i), :);
end

figure(10);
plot(node_density, max_d, 'ko');
title("Maximum displacement");
xlabel("Node density (#/m)");
ylabel("Displacement magnitude (m)");

figure(11);
plot(node_density, max_s_v, 'ko', node_density, max_s_t, 'kx');
title("Maximum stress");
xlabel("Node density (#/m)");
ylabel("Maximum stress (Pa)");


% make boundary conditions for finer mesh
function res = boundary_cond(node_density, w)
    
    node_density = round(node_density);
    x_num = 20 * node_density + 1;
    y_num = 1 * node_density + 1;
    nodes = zeros(x_num * y_num, 2);
    for i = 1:x_num
        for j = 1:y_num
            nodes((i-1)*y_num+j, :) = [ ...
                (i-1)/(x_num-1)*20, (j-1)/(y_num-1)];
        end
    end
    element = delaunay(nodes);
    dl = 1 / node_density;
    df = dl * (-w);
    bcond = [1:2*x_num*y_num; zeros(1, 2*x_num*y_num)];
    for i = y_num:y_num:size(nodes, 1)
        bcond(2, 2*i) = df;
    end
    bcond = bcond(:, 2*y_num+1:2*(x_num-1)*y_num);
    res = {element, nodes, bcond};

end

% obtain stress for each node
function stress = node_stress(nodes, elements, element_stress)
    
    numbers = zeros(size(nodes, 1), 1);
    stress = zeros(size(nodes, 1), 3);
    for i = 1:size(elements, 1)
        node_index = elements(i, :);
        numbers(node_index, 1) = numbers(node_index, 1) + [1; 1; 1];
        stress(node_index, :) = ...
            stress(node_index, :) + element_stress(i, :);
    end
    stress(:, 1) = stress(:, 1) ./ numbers;
    stress(:, 2) = stress(:, 2) ./ numbers;
    stress(:, 3) = stress(:, 3) ./ numbers;
        
end

% obtain von mises stress / plane strain
function stress_von_mises = von_mises(stress, nu)

    s_x = stress(:, 1);
    s_y = stress(:, 2);
    s_z = nu * (s_x + s_y);
    t_xy = stress(:, 3);
    stress_von_mises = sqrt( ...
        0.5 * (s_y - s_z).^2 ...
        + 0.5 * (s_z - s_x).^2 ...
        + 0.5 * (s_x - s_y).^2 ...
        + 3 * t_xy.^2);

end

% obtain tresca stress / plane strain
function stress_tresca = tresca(stress, nu)

    stress_tresca = zeros(size(stress, 1), 1);
    s_x = stress(:, 1);
    s_y = stress(:, 2);
    t_xy = stress(:, 3);
    s_z = nu * (s_x + s_y);
    for i = 1:size(stress, 1)
        s = [s_x(i, 1), t_xy(i, 1), 0;
            t_xy(i, 1), s_y(i, 1),  0;
            0,          0,          s_z(i, 1)];
        s_p = eig(s);
        stress_tresca(i, 1) = max(abs(s_p - s_p([2, 3, 1], 1)));
    end

end

% obtain stress for each element
function stress = element_stress(nodes, elements, dx, dy, D)

    stress = zeros(size(elements, 1), 3);

    for i = 1:size(elements, 1)
        
        element_nodes = nodes(elements(i, :), :);
        area = polyarea(element_nodes(:, 1), element_nodes(:, 2));
        b1 = element_nodes(2, 2) - element_nodes(3, 2);
        b2 = element_nodes(3, 2) - element_nodes(1, 2);
        b3 = element_nodes(1, 2) - element_nodes(2, 2);
        c1 = element_nodes(3, 1) - element_nodes(2, 1);
        c2 = element_nodes(1, 1) - element_nodes(3, 1);
        c3 = element_nodes(2, 1) - element_nodes(1, 1);
        B = [
            b1, 0 , b2, 0 , b3, 0 ;
            0 , c1, 0 , c2, 0 , c3;
            c1, b1, c2, b2, c3, b3
        ] / (2 * area);
        d = [dx(elements(i, :), 1), dy(elements(i, :), 1)];
        d = reshape(d', [], 1);
        stress(i, :) = D * B * d;

    end

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
 
    % iteration for each element
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

% solve for K_global, boundary condition
function [d, f] = solve_for(K_global, boundary_cond)

    % Reduced stiffness matrix
    K_reduced = sparse(K_global(boundary_cond(1, :), boundary_cond(1, :)));
    
    % deflection and force at each node
    d = zeros(size(K_global, 1), 1);
    f = zeros(size(K_global, 1), 1);
    
    % reduced f
    f(boundary_cond(1, :), 1) = boundary_cond(2, :)';
    
    % obtain d from solving the matrix equation
    d(boundary_cond(1, :), 1) = K_reduced \ boundary_cond(2, :)';
    
end
