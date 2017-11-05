boundary = [
    0, -2;
    0, 2;
    -5, 1.5;
    -10, 1;
    -10, -1;
    -5, -1.5
    ];

E = 200*10^9;
thickness = 1;

mesh = generate_mesh(boundary, 0.3);

elements = delaunay(mesh);

K_g = zeros(2 * size(mesh, 1));

element_num = size(elements, 1);
for i = element_num:-1:1
    element = elements(i, :);
    area = polyarea(mesh(element, 1), mesh(element, 2));
    if area < 10^-5
        elements(i, :) = [];
    end
end
element_num = size(elements, 1);

for i = 1:element_num
    element = elements(i, :);
    area = polyarea(mesh(element, 1), mesh(element, 2));
    index = [ ...
        2 * elements(i, 1) - 1, ...
        2 * elements(i, 1), ...
        2 * elements(i, 2) - 1, ...
        2 * elements(i, 2), ...
        2 * elements(i, 3) - 1, ...
        2 * elements(i, 3) ...
        ];
    b1 = mesh(element(1, 2), 2) - mesh(element(1, 3), 2);
    b2 = mesh(element(1, 3), 2) - mesh(element(1, 1), 2);
    b3 = mesh(element(1, 1), 2) - mesh(element(1, 2), 2);
    c1 = mesh(element(1, 3), 1) - mesh(element(1, 2), 1);
    c2 = mesh(element(1, 1), 1) - mesh(element(1, 3), 1);
    c3 = mesh(element(1, 2), 1) - mesh(element(1, 1), 1);
    B = [
        b1, 0, b2, 0, b3, 0;
        0, c1, 0, c2, 0, c3;
        c1, b1, c2, b2, c3, b3
        ] / (2 * area);
    D = [
        1, 0.3, 0;
        0.3, 1, 0;
        0, 0, 0.35
        ] * E / 0.91;
    K_element = thickness * area * B' * D * B;
    K_g(index, index) = K_g(index, index) + K_element;
end

bcond = generate_bcond(mesh);

d = zeros(size(K_g, 1), 1);

K_reduced = K_g(bcond(1, :), bcond(1, :));
f_reduced = bcond(2, :)';
d(bcond(1, :), 1) = K_reduced \ f_reduced;

dx = d(1:2:size(d, 1), 1);
dy = d(2:2:size(d, 1), 1);

triplot(elements, mesh(:, 1) + 20000 * dx, mesh(:, 2) + 20000 * dy);


axis equal

function bcond = generate_bcond(node)
    P1 = 100000;
    P2 = 200000;
    bcond = [];
    num = size(node, 1);
    for i = 1:num
        if node(i, 1) ~= 0
            if node(i, :) == [-5, 1.5]
                bcond = [bcond, [2*i-1, 2*i; 0, -P2]];
            elseif node(i, :) == [-10, 1]
                bcond = [bcond, [2*i-1, 2*i; 0, -P1]];
            else
                bcond = [bcond, [2*i-1, 2*i; 0, 0]];
            end
        end
    end
end

function xy = generate_mesh(boundary, len)
    xy = boundary;
    xsize = 10;
    xnum = ceil(5 / len) * 2;
    for i = 0:xnum
        xx = -10 * i / xnum;
        ylen = 4 + 0.2 * xx;
        ynum = ceil(ylen / len);
        for j = 0:ynum
            yup = 0.5 * ylen;
            ydown = -0.5 * ylen;
            yy = (j * yup + (ynum - j) * ydown) / ynum;
            
            xy = [xy; xx, yy];
        end
    end
    xy = unique(xy, "row");
end