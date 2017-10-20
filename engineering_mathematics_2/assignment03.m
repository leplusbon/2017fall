% Engineering Mathematics 2
% Assignment 3
% 2014-10773 ¿Â¿Œ±‘

% problem 2

dx = 0.01;      % dx for plotting
nsize = 20;     % # of terms of the fourier series

% x axis
x = -5:dx:5;

% coefficients
% b_i's are all zero
a0 = pi / 2;
a = zeros(nsize, 1);
for i = 1:2:nsize
    a(i) = 4 / i / i / pi;
end

fourier(1, :) = a0 * ones(size(x));
for i = 1:nsize
    fourier(i + 1, :) = a(i) * cos(i * x);
end

% given function f in [-5, 5]
f_given = abs(abs(x) - pi);

% partial sum up to 5th, 20th terms
f_5 = sum(fourier(1:6, :), 1);
f_20 = sum(fourier(1:21, :), 1);

% plot results
plot(x, f_given, x, f_5, x, f_20);

