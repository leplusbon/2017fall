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
    1, 2, 7, 8;
    2, 3, 8, 7;
    3, 4, 9, 8;
    4, 5, 10, 9;
    6, 7, 12, 11;
    7, 8, 13, 12;
    8, 9, 14, 13;
    9, 10, 15, 14
    ];
    
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
        
        K_g(K_index, K_index) = ...
            K_g(K_index, K_index) + get_k_element(element, E, v, t);
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

end

function [mesh, elements] = generate_mesh(boundary, mesh_size)
    n = round(2 / mesh_size);
    x_coord = 0:1:(n*4);
    y_coord = (-n*2):1:(n*2);
    
    x_coord = x_coord / n;
    y_coord = y_coord / n;
    
    
end

function [max_absd, max_absd_node, max_stressv, max_stressv_node] ...
    = postprocess(dx, dy, mesh, elements)

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
    ttss = [
        -alpha, -alpha;
        +alpha, -alpha;
        +alpha, +alpha;
        -alpha, +alpha
        ];

    for i = 1:4
        tt = ttss(i, 1);
        ss = ttss(i, 2);
        
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
        k_element = k_element + B' * D * B * J * t * 1 * 1;
        
    end
    
end