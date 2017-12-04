%% given values

% thickness (m)
t = 100;

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

elem = delaunay(initial_node);
triplot(elem, initial_node(:, 1), initial_node(:, 2), 'k');
axis equal

%%
