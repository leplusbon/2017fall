

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

    %stress = element_stress(nodes, elements, dx, dy, D);
    stress = node_stress(nodes, elements, dx, dy, D);
    s_v = von_mises(stress, nu);
    s_t = tresca(stress, nu);
    [max_d(i), node_num_d(i)] = max(sqrt(dx.^2 + dy.^2));
    [max_s_v(i), node_num_s_v(i)] = max(s_v);
    [max_s_t(i), node_num_s_t(i)] = max(s_t);
    loc_d(i, :) = nodes(node_num_d(i), :);
    loc_s_v(i, :) = nodes(node_num_s_v(i), :);
    loc_s_t(i, :) = nodes(node_num_s_t(i), :);
end
% 
% xx = nodes(:, 1) + 100000 * dx;
% yy = nodes(:, 2) + 100000 * dy;
% 
% xx = flipud(reshape(xx, sqrt(size(elements, 1) / 20) + 1, ...
%     sqrt(size(elements, 1) / 20) * 20 + 1));
% yy = flipud(reshape(yy, sqrt(size(elements, 1) / 20) + 1, ...
%     sqrt(size(elements, 1) / 20) * 20 + 1));
% s_v = flipud(reshape(s_v, sqrt(size(elements, 1) / 20) + 1, ...
%     sqrt(size(elements, 1) / 20) * 20 + 1));
% figure(9);
% pcolor(xx, yy, s_v);
% colormap(jet);
% title("Deformed shape, von Mises stress (Pa) - 8000 square elements");
% xlabel("x (m)");
% ylabel("y (m)");
% axis equal
% 
% figure(10);
% plot(node_density, max_d, 'ko');
% title("Maximum displacement");
% xlabel("Node density (#/m)");
% ylabel("Displacement magnitude (m)");
% 
% figure(11);
% plot(node_density, max_s_v, 'ko', node_density, max_s_t, 'kx');
% title("Maximum stress");
% xlabel("Node density (#/m)");
% ylabel("Maximum stress (Pa)");

% make boundary conditions for finer square mesh
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
    elements = zeros((x_num-1) * (y_num-1), 4);
    index = 1;
    for i = 1:x_num-1
        for j = 1:y_num-1
            elements(index, :) = [ ...
                (i - 1) * y_num + j, ...
                i * y_num + j, ...
                i * y_num + j + 1, ...
                (i - 1) * y_num + j + 1
                ];
            index = index + 1;
        end
    end
    dl = 1 / node_density;
    df = dl * (-w);
    bcond = [1:2*x_num*y_num; zeros(1, 2*x_num*y_num)];
    for i = y_num:y_num:size(nodes, 1)
        bcond(2, 2*i) = df;
    end
    bcond = bcond(:, 2*y_num+1:2*(x_num-1)*y_num);
    res = {elements, nodes, bcond};

end

% obtain stress for each node
function stress = node_stress(nodes, elements, dx, dy, D)
    
    num = zeros(size(nodes, 1), 1);
    stress = zeros(size(nodes, 1), 3);
    
    for elm_index = 1:size(elements, 1)
        
        x = nodes(elements(elm_index, :), 1);
        y = nodes(elements(elm_index, :), 2);
        
        alpha = 1;

        % perform gaussian quadrature(4 points) to get K_element
        sstt = [
            -alpha, -alpha;
            +alpha, -alpha;
            +alpha, +alpha;
            -alpha, +alpha
            ];

        for i = 1:4
            ss = sstt(i, 1);
            tt = sstt(i, 2);

            % determine |J|
            J = 0.125 * x' * [
                0,      1-tt,   tt-ss,  ss-1;
                tt-1,   0,      ss+1,   -ss-tt;
                ss-tt,  -ss-1,  0,      tt+1;
                1-ss,   ss+tt,  -tt-1,  0
                ] * y;

            % determine B
            a = 0.25 * [ss-1, -1-ss, 1+ss, 1-ss] * y;
            b = 0.25 * [tt-1, 1-tt, 1+tt, -1-tt] * y;
            c = 0.25 * [tt-1, 1-tt, 1+tt, -1-tt] * x;
            d = 0.25 * [ss-1, -1-ss, 1+ss, 1-ss] * x;

            Ns(1) = 0.25 * (tt - 1);
            Nt(1) = 0.25 * (ss - 1);
            Ns(2) = 0.25 * (1 - tt);
            Nt(2) = 0.25 * (-1 - ss);
            Ns(3) = 0.25 * (1 + tt);
            Nt(3) = 0.25 * (1 + ss);
            Ns(4) = 0.25 * (-1 - tt);
            Nt(4) = 0.25 * (1 - ss);

            B = zeros(3, 8);
            for j = 1:4
                B(:, (2*j-1):(2*j)) = [
                    a*Ns(j)-b*Nt(j),    0;
                    0,                  c*Nt(j)-d*Ns(j);
                    c*Nt(j)-d*Ns(j),    a*Ns(j)-b*Nt(j)
                    ];
            end
            B = B / J;
            
            node_index = elements(elm_index, i);
            d = [dx(elements(elm_index, :), 1), dy(elements(elm_index, :), 1)];
            d = reshape(d', [], 1);
            stress(node_index, :) = stress(node_index, :) + (D * B * d)';
            num(node_index, 1) = num(node_index, 1) + 1;
            
        end
    end
    
    for i = 1:3
        stress(:, i) = stress(:, i) ./ num;
    end
    
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

% obtain element stiffness matrix for rectangular element
function K_element = element_stiffness(nodes, thickness, D)
 
    K_element = zeros(8, 8);

    x = nodes(:, 1);
    y = nodes(:, 2);

    % gaussian points for integral
    alpha = 1 / sqrt(3);
    
    % perform gaussian quadrature(4 points) to get K_element
    sstt = [
        -alpha, -alpha;
        +alpha, -alpha;
        +alpha, +alpha;
        -alpha, +alpha
        ];

    for i = 1:4
        ss = sstt(i, 1);
        tt = sstt(i, 2);
        
        % determine |J|
        J = 0.125 * x' * [
            0,      1-tt,   tt-ss,  ss-1;
            tt-1,   0,      ss+1,   -ss-tt;
            ss-tt,  -ss-1,  0,      tt+1;
            1-ss,   ss+tt,  -tt-1,  0
            ] * y;
        
        % determine B
        a = 0.25 * [ss-1, -1-ss, 1+ss, 1-ss] * y;
        b = 0.25 * [tt-1, 1-tt, 1+tt, -1-tt] * y;
        c = 0.25 * [tt-1, 1-tt, 1+tt, -1-tt] * x;
        d = 0.25 * [ss-1, -1-ss, 1+ss, 1-ss] * x;
        
        Ns(1) = 0.25 * (tt - 1);
        Nt(1) = 0.25 * (ss - 1);
        Ns(2) = 0.25 * (1 - tt);
        Nt(2) = 0.25 * (-1 - ss);
        Ns(3) = 0.25 * (1 + tt);
        Nt(3) = 0.25 * (1 + ss);
        Ns(4) = 0.25 * (-1 - tt);
        Nt(4) = 0.25 * (1 - ss);
        
        B = zeros(3, 8);
        for j = 1:4
            B(:, (2*j-1):(2*j)) = [
                a*Ns(j)-b*Nt(j),    0;
                0,                  c*Nt(j)-d*Ns(j);
                c*Nt(j)-d*Ns(j),    a*Ns(j)-b*Nt(j)
                ];
        end
        B = B / J;
        
        % perform guassian quadrature
        K_element = K_element + J * thickness * (B' * D * B);
    end
    
end
 
% obtain global stiffness matrix for rectangular mesh
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
