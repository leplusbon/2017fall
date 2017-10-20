%% 4

a = 174.56;
b = 0.001381;
R = 293.8;
p = 2000000;
T = 300;

v_0 = -1;
v_1 = R * T / p;

while v_1 ~= v_0
    v_0 = v_1;
    y_0 = p * v_0^3 - (p * b + R * T) * v_0^2 + a * v_0 - a * b;
    v_p = 3 * p * v_0^2 - 2 * (p * b + R * T) * v_0 + a;
    v_1 = v_0 - y_0 / v_p;
end

disp("[5] v (m3)");
disp(v_1);

%% 5

f_0 = -1;
f_1 = 0.01;

while f_1 ~= f_0
    f_0 = f_1;
    f_1 = 1 / (11.2 * log10(f_0))^2;
end

disp("[6] f");
disp(f_1);

%% 6
