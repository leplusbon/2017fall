
boundary_initial = [
    0, -2;
    0, 2;
    -5, 1.5;
    -10, 1;
    -10, -1;
    -5, -1.5
    ];

tri = delaunay(boundary_initial);

xy_initial = generate_initial_boundary(boundary_initial, 0.5);


el = generate_mesh(boundary_initial, xy_initial, delaunay(xy_initial), 100, cosd(18));

function xy = generate_initial_boundary(boundary_initial, dist)
    boundary_size = size(boundary_initial, 1);
    xy = boundary_initial;
    for i = 1:boundary_size
        length = norm(boundary_initial(i, :) - boundary_initial(mod(i, boundary_size) + 1, :));
        node_num = round(length / dist);
        xy_addition = zeros(node_num - 1, 2);
        for j = 1:(node_num - 1)
            xy_addition(j, :) = ...
                (j * boundary_initial(i, :) + ...
                (node_num - j) * boundary_initial(mod(i, boundary_size) + 1, :)) / ...
                node_num;
        end
        xy = [xy; xy_addition];
    end
end

function elements = generate_mesh(boundary_initial, xy_initial, tri_initial, max_area, max_cosine)
    
    xy = insert_vertices(xy_initial, tri_initial, max_area, max_cosine, boundary_initial);
    tri = delaunay(xy(:, 1), xy(:, 2));

    for i = 1:100
        disp(i);
        num = size(xy, 1);
        xy = insert_vertices(xy, tri, max_area, max_cosine, boundary_initial);
        if num == size(xy, 1)
            break;
        end
        tri = delaunay(xy(:, 1), xy(:, 2));
    end
    if i > 99
        disp("failure generating mesh");
        elements = xy;
        return;
    end
    elements = xy;
end

function poor = is_poor(x, y, max_area, max_cosine)
    poor = 0;
    if polyarea(x, y) > max_area
        poor = 1;
        return;
    end
    a2 = (x(2, 1) - x(1, 1))^2 + (y(2, 1) - y(1, 1))^2;
    b2 = (x(3, 1) - x(2, 1))^2 + (y(3, 1) - y(2, 1))^2;
    c2 = (x(1, 1) - x(3, 1))^2 + (y(1, 1) - y(3, 1))^2;
    if a2 + b2 - c2 > 2 * sqrt(a2 * b2) * max_cosine
        poor = 1;
        return;
    elseif b2 + c2 - a2 > 2 * sqrt(b2 * c2) * max_cosine
        poor = 1;
        return;
    elseif c2 + a2 - b2 > 2 * sqrt(c2 * a2) * max_cosine
        poor = 1;
        return;
    end
end

function inside = is_inside_boundary(x, y, xy_initial)
    inside = inpolygon(x, y, xy_initial(:, 1), xy_initial(:, 2));
end

function encroached = is_encroached(xy, cc)
    a2 = (xy(1, 1) - xy(2, 1))^2 + (xy(1, 2) - xy(2, 2))^2;
    b2 = (xy(1, 1) - cc(1, 1))^2 + (xy(1, 2) - cc(1, 2))^2;
    c2 = (xy(2, 1) - cc(1, 1))^2 + (xy(2, 2) - cc(1, 2))^2;
    
    if b2 + c2 - a2 < 0
        encroached = 1;
    else
        encroached = 0;
    end 
end

function xy_new = insert_vertices(xy, tri, max_area, max_cosine, xy_initial)
    element_num = size(tri, 1);
    for i = 1:element_num
        xx = xy(tri(i, :), 1);
        yy = xy(tri(i, :), 2);
        if is_poor(xx, yy, max_area, max_cosine) == 1
            tr = triangulation([1, 2, 3], xx, yy);
            cc = circumcenter(tr);
            if is_inside_boundary(cc(1, 1), cc(1, 2), xy_initial)
                if is_encroached([xx([1, 2], 1), yy([1, 2], 1)], cc) == 1
                    xy = [xy; 0.5 * (xx(2, 1) + xx(1, 1)), ...
                        0.5 * (yy(2, 1) + yy(1, 1))];
                elseif is_encroached([xx([2, 3], 1), yy([2, 3], 1)], cc) == 1
                    xy = [xy; 0.5 * (xx(3, 1) + xx(2, 1)), ...
                        0.5 * (yy(3, 1) + yy(2, 1))];
                elseif is_encroached([xx([3, 1], 1), yy([3, 1], 1)], cc) == 1
                    xy = [xy; 0.5 * (xx(1, 1) + xx(3, 1)), ...
                        0.5 * (yy(1, 1) + yy(3, 1))];
                else
                    xy = [xy; cc];
                end
            end
        end
    end
    xy_new = unique(xy, 'rows');
end