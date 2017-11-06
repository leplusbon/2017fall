initial_mesh = [
    0, -1;
    1, -1;
    2, -1;
    3, -1;
    4, -1;
    0, 0;
    1, 0;
    2, 0;
    3, 0;
    4, 0;
    0, 1;
    1, 1;
    2, 1;
    3, 1;
    4, 1
    ];

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
    

function [dx, dy] = solve(mesh, elements, bcond)

end

function bcond = boundary_cond(mesh)

end

function mesh = generate_mesh(boundary, mesh_size)

end

function [max_absd, max_absd_node, max_stressv, max_stressv_node] ...
    = postprocess(dx, dy, mesh, elements)

end

