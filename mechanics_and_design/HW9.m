% given nodes
initial_mesh = [
    0, -2;
    2, -2;
    4, -2;
    6, -2;
    8, -2;
    0, 0;
    2, 0;
    4, 0;
    6, 0;
    8, 0;
    0, 2;
    2, 2;
    4, 2;
    6, 2;
    8, 2
    ];

% given elements
initial_elements = [
    1, 2, 7, 6;
    2, 3, 8, 7;
    3, 4, 9, 8;
    4, 5, 10, 9;
    6, 7, 12, 11;
    7, 8, 13, 12;
    8, 9, 14, 13;
    9, 10, 15, 14
    ];

[dx, dy] = solve(initial_mesh, initial_elements, boundary_cond(initial_mesh));

max_stress_node = 0;
max_stress_v = 0;
max_absd_node = 0;
max_absd = 0;
for i = 1:15
    if max_stress_v < von_mises(get_stress(initial_mesh, initial_elements, i, dx, dy, 200*10^9, 0.3))
        max_stress_node = i;
        max_stress_v = von_mises(get_stress(initial_mesh, initial_elements, i, dx, dy, 200*10^9, 0.3));
    end
    if max_absd < sqrt(dx(i, 1)^2 + dy(i, 1)^2)
        max_absd_node = i;
        max_absd = sqrt(dx(i, 1)^2 + dy(i, 1)^2);
    end
end

disp("(a)");
disp("point of maximum deflection (m) :");
disp(initial_mesh(max_absd_node, :));
disp("maximum deflection (m) :");
disp(max_absd);

disp("(b)");
disp("point of maximum von mises stress (m) :");
disp(initial_mesh(max_stress_node, :));
disp("maximum von mises stress (Pa) :");
disp(max_stress_v);

max_stress_v = [];
num_element = [];
max_absd = [];

for i = 2:2:40
    mesh_size = 4 / i;
    fprintf("solving for mesh_size = %f\n", mesh_size);
    [mesh, elements] = generate_mesh(mesh_size);
    bcond = boundary_cond(mesh);
    [dx, dy] = solve(mesh, elements, bcond);

    stress_v = zeros(size(mesh, 1), 1);
    for j = 1:size(mesh, 1)
        stress_v(j, 1) = von_mises(get_stress(mesh, elements, j, dx, dy, 200*10^9, 0.3));
    end

    xx = mesh(:, 1) + 30000 * dx;
    yy = mesh(:, 2) + 30000 * dy;

    xx = reshape(xx, sqrt(2 * size(elements, 1)) + 1, sqrt(0.5 * size(elements, 1)) + 1)';
    yy = reshape(yy, sqrt(2 * size(elements, 1)) + 1, sqrt(0.5 * size(elements, 1)) + 1)';
    stress_v = reshape(stress_v, sqrt(2 * size(elements, 1)) + 1, sqrt(0.5 * size(elements, 1)) + 1)';
    max_stress_v = [max_stress_v, max(stress_v(:))];
    max_absd = [max_absd, sqrt(max((dx .* dx + dy .* dy)))];
    num_element = [num_element, size(elements, 1)];
    
end

figure(1);
plot(num_element, max_absd);
title("(c) maximum deflection");
xlabel("# of elements");
ylabel("maximum deflection (m)");

figure(2);
plot(num_element, max_stress_v);
title("(d) maximum von mises stress")
xlabel("# of elements");
ylabel("maximum von mises stress (Pa)");

figure(3);
pcolor(xx, yy, stress_v);
colormap(jet);
title("(e) deformed shape(3200 elements, deflection exaggerated, x30000)");
xlabel("x (m)");
ylabel("y (m)");
axis equal

function stress_v = von_mises(stress)
    sigma1 = stress(1, 1);
    sigma2 = stress(2, 1);
    tau = stress(3, 1);
    
    stress_v = sqrt(0.5 * ( ...
        (sigma1 - sigma2)^2 + ...
        sigma1^2 + ...
        sigma2^2 + ...
        6 * tau^2 ...
    ));
    
end

% solve for the system
function [dx, dy] = solve(mesh, elements, bcond)
    
    % material properties
    E = 200 * 10^9;
    v = 0.3;
    t = 0.5;

    node_num = size(mesh, 1);
    K_g = zeros(2 * node_num);
    
    element_num = size(elements, 1);
    for i = 1:element_num
        element_index = elements(i, :);
        element = mesh(element_index, :);
        K_index = [2*element_index-1; 2*element_index];
        K_index = (K_index(:))';
        k_elem = get_k_element(element, E, v, t);
        K_g(K_index, K_index) = ...
            K_g(K_index, K_index) + k_elem;
    end
    
    d = zeros(2 * node_num, 1);
    
    % solve by performing d = K^(-1) * f
    K_reduced = K_g(bcond(1, :), bcond(1, :));
    f_reduced = bcond(2, :)';
    d(bcond(1, :), 1) = K_reduced \ f_reduced;

    dx = d(1:2:size(d, 1), 1);
    dy = d(2:2:size(d, 1), 1);
    
end

function bcond = boundary_cond(mesh)
    node_num = size(mesh, 1);
    bcond = [1:(node_num*2); zeros(1, node_num*2)];
    
    nl = (-6 + sqrt(8 * node_num + 1)) / 8;
    P = -8 * 10000 * 0.5 / nl;
    for i = 1:node_num
        if mesh(i, 2) == 2
            bcond(2, 2 * i) = P;
        end
    end
    for i = node_num:(-1):1
        if mesh(i, 1) == 0
            bcond(:, 2 * i) = [];
            bcond(:, 2 * i - 1) = [];
        end
    end
    
end

function [mesh, elements] = generate_mesh(mesh_size)

    n = round(2 / mesh_size);
    x_coord = 0:1:(n*4);
    y_coord = (-n):1:(n);
    
    x_coord = 2 * x_coord / n;
    y_coord = 2 * y_coord / n;
    
    node_num = (n * 4 + 1) * (n * 2 + 1);
    mesh = zeros(node_num, 2);
    
    for i = 1:(n*4+1)
        for j = 1:(n*2+1)
            mesh((j - 1) * (n * 4 + 1) + i, :) = [
                x_coord(i), y_coord(j)
                ];
        end
    end
    
    elements = zeros(8 * n * n, 4);
    for i = 1:(n*4)
        for j = 1:(n*2)
            elements((j - 1) * n * 4 + i, :) = [ ...
                (j - 1) * (n * 4 + 1) + i, ...
                (j - 1) * (n * 4 + 1) + i + 1, ...
                j * (n * 4 + 1) + i + 1, ...
                j * (n * 4 + 1) + i  ...
                ];
        end
    end
    
end

% get element stiffness matrix for 4-node element
function k_element = get_k_element(element, E, v, t)
    
    k_element = zeros(8, 8);

    x = element(:, 1);
    y = element(:, 2);
    
    % D matrix for plane stress condition
    D = [
        1, v, 0;
        v, 1, 0;
        0, 0, (1-v)/2
        ] * E / (1 - v^2);
    
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
        k_element = k_element + J * t * (B' * D * B);
    end
    
end

function stress = get_stress(mesh, elements, node_index, dx, dy, E, v)

    [row_index, column_index] = find(elements == node_index);
    n = size(row_index, 1);
    
    D = [
        1, v, 0;
        v, 1, 0;
        0, 0, (1-v)/2
        ] * E / (1 - v^2);
    
%     sstt = [
%         -1, -1;
%         +1, -1;
%         +1, +1;
%         -1, +1
%         ];

    alpha = 0;
    
    sstt = [
        -alpha, -alpha;
        +alpha, -alpha;
        +alpha, +alpha;
        -alpha, +alpha
        ];
    
    stress = zeros(3, 1);
    
    for i = 1:n
        
        element_dx = dx(elements(row_index(i, 1), :), 1)';
        element_dy = dy(elements(row_index(i, 1), :), 1)';
        
        d_element = [element_dx; element_dy];
        
        d_element = d_element(:);
        
        x = mesh(elements(row_index(i, 1), :), 1);
        y = mesh(elements(row_index(i, 1), :), 2);
        
        ss = sstt(column_index(i, 1), 1);
        tt = sstt(column_index(i, 1), 2);
        
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
        
        stress = stress + D * B * d_element;
        
    end
    
    stress = stress / n;
    
end

