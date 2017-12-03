%%
for i = 1:21
    boundary(i, :) = [i-1, 0];
    boundary(43-i, :) = [i-1, 1];
end

elem = delaunay(boundary);
triplot(elem, boundary(:, 1), boundary(:, 2));
axis equal
%%
initial_elements = [
    1, 2, 6;
    2, 3, 6;
    3, 5, 6;
    3, 4, 5
    ];

mesh = boundary;

elements = initial_elements;

[dx, dy] = solve(mesh, elements);
[max_absd, maxdelem, max_stressv, maxstresselem] ...
    = postprocess(dx, dy, mesh, elements);

disp("(a)");
disp("location of maximum deflection (m) : ");
disp(mesh(maxdelem, :));
disp("maximum deflection (m) : ");
disp([dx(maxdelem, 1), dy(maxdelem, 1)]);

disp("(b)");
disp("location of maximum von mises stress (m) : ");
disp(mesh(maxstresselem, :));
disp("maximum von mises stress (Pa) : ");
disp(max_stressv);


result = [4, max_absd, max_stressv];
for i = 4:2:36
    dl = 5 / i;
    mesh = generate_mesh(boundary, dl);
    elements = make_elements(mesh);
    numelem = size(elements, 1);
    fprintf("running fem for %d elements\n", numelem);
    [dx, dy] = solve(mesh, elements);
    [max_absd, maxdelem, max_stressv, maxstresselem] ...
        = postprocess(dx, dy, mesh, elements);
    result = [result; [numelem, max_absd, max_stressv]];
end

figure(1);
plot(result(:, 1), result(:, 2));
title("(c) maximum deflection");
xlabel("# of elements");
ylabel("maximum deflection (m)");

figure(2);
plot(result(:, 1), result(:, 3));
title("(d) maximum von mises stress");
xlabel("# of elements");
ylabel("maximum von mises stress (Pa)");

figure(3);
triplot(elements, mesh(:, 1) + 15000 * dx, mesh(:, 2) + 15000 * dy);
title("(e) deflection shape (deflection exaggerated, x15000, 3168 elements)");
xlabel("x (m)");
ylabel("y (m)");
axis equal

function elements = make_elements(mesh)
    elements = delaunay(mesh);
    element_num = size(elements, 1);
    for i = element_num:-1:1
        element = elements(i, :);
        area = polyarea(mesh(element, 1), mesh(element, 2));
        if area < 10^-5
            elements(i, :) = [];
        end
    end
end

function [dx, dy] = solve(mesh, elements)
    E = 200*10^9;
    thickness = 1;
    
    K_g = zeros(2 * size(mesh, 1));

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

end

function [max_absd, maxdelem, max_stressv, maxstresselem] ...
    = postprocess(dx, dy, mesh, elements)

    E = 200*10^9;
    max_absd = 0;
    maxdelem = 0;
    max_stressv = 0;
    
    d = [dx, dy]';
    d = d(:);
    
    for i = 1:size(dx, 1)
        absd = sqrt(dx(i, 1)^2 + dy(i, 1)^2);
        if max_absd < absd
            maxdelem = i;
            max_absd = absd;
        end
        [elm, temp] = find(elements == i);
        elem_on_node = size(elm, 1);
        stress = zeros(3, 1);
        for j = 1:elem_on_node
            element = elements(elm(j, 1), :);
            area = polyarea(mesh(element, 1), mesh(element, 2));
            index = [ ...
                2 * element(1, 1) - 1, ...
                2 * element(1, 1), ...
                2 * element(1, 2) - 1, ...
                2 * element(1, 2), ...
                2 * element(1, 3) - 1, ...
                2 * element(1, 3) ...
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
            stress = stress + D * B * d(index, 1);
            

        end
        stress = stress / elem_on_node;
        stress_v = sqrt( ...
            ((stress(1, 1) - stress(2, 1))^2 + ...
            stress(1, 1)^2 + ...
            stress(2, 1)^2 + ...
            6 * stress(3, 1)^2) / 2);

        if max_stressv < stress_v
            max_stressv = stress_v;
            maxstresselem = i;
        end
    end
end

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
        xx = -xsize * i / xnum;
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

