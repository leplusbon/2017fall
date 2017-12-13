Estar = 219780219780;
x = 0:0.01:20;
dy = -100000/(24*Estar*8.33333) * (x .* x) .* ((20-x) .* (20-x));
tresca = 100000*0.25/12/8.333333*abs(400-120*x+6*x.^2);
vonmises = tresca * sqrt(1-0.3+0.09);

tresca_fem_x = 0.05:0.1:19.95;
vonmises_fem_x = 0.05:0.1:19.95;
tresca_fem_y = zeros(200, 1);
vonmises_fem_y = zeros(200, 1);
for i=1:200
    index = (10*i-9):(10*i);
    tresca_fem_y(i) = mean(tresca_fem(index, 1));
    vonmises_fem_y(i) = mean(vonmises_fem(index, 1));
end

x_tri = zeros(601, 1);
y_tri = zeros(601, 1);
x_rec = zeros(601, 1);
y_rec = zeros(601, 1);
vonmises_tri = zeros(601, 1);
tresca_tri = zeros(601, 1);
vonmises_rec = zeros(601, 1);
tresca_rec = zeros(601, 1);
for i=1:601
    x_tri(i) = (i-1) * 1/30;
    x_rec(i) = (i-1) * 1/30;
    for j = 1:31
        vonmises_tri(i) = vonmises_tri(i) + s_v_tri((i-1)*31+j, 1);
        vonmises_rec(i) = vonmises_rec(i) + s_v_rec((i-1)*31+j, 1);
        tresca_tri(i) = tresca_tri(i) + s_t_tri((i-1)*31+j, 1);
        tresca_rec(i) = tresca_rec(i) + s_t_rec((i-1)*31+j, 1);
        y_tri(i) = y_tri(i) + dy_tri((i-1)*31+j, 1);
        y_rec(i) = y_rec(i) + dy_rec((i-1)*31+j, 1);
    end
end
y_tri = y_tri / 31;
y_rec = y_rec / 31;
vonmises_tri = vonmises_tri / 31;
vonmises_rec = vonmises_rec / 31;
tresca_tri = tresca_tri / 31;
tresca_rec = tresca_rec / 31;

figure(1);
plot(x, dy, 'k-', x_fem(:, 1), x_fem(:, 2), 'k--', x_tri, y_tri, 'k:', x_rec, y_rec, 'k-.');
legend('Analytic', 'FEM(HyperWorks)', 'FEM(triangular)', 'FEM(rectangular)');
xlabel("x (m)");
ylabel("y (m)");
title("Deflection in y direction along the neutral axis");

figure(2);
plot(x, vonmises, 'k-', vonmises_fem_x, vonmises_fem_y, 'k--', x_tri, vonmises_tri, 'k:', x_rec, vonmises_rec, 'k-.');
legend('Analytic', 'FEM(HyperWorks)', 'FEM(triangular)', 'FEM(rectangular)');
xlabel("x (m)");
ylabel("\sigma_v (Pa)");
title("Von Mises stress");

figure(3);
plot(x, tresca, 'k-', tresca_fem_x, tresca_fem_y, 'k--', x_tri, tresca_tri, 'k:', x_rec, tresca_rec, 'k-.');
legend('Analytic', 'FEM(HyperWorks)', 'FEM(triangular)', 'FEM(rectangular)');
xlabel("x (m)");
ylabel("\sigma_{Tresca} (Pa)");
title("Tresca stress");