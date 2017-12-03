table2 = [
    331, 0, 43000;
    276, 0, 112000;
    241, 0, 172000;
    207, 0, 231000;
    190, 0, 546000;
    179, 0, 1165000
    ];

table3 = [
    293, 592, 45000;
    241, 646, 90000;
    207, 668, 160000;
    174, 685, 700000;
    778, -130, 4500;
    594, -312, 132000;
    529, -354, 540000
    ];

s_u = 1233;
s_fB = 1717;

data = [table2; table3];

s_a = data(:, 1);
s_m = data(:, 2);
N_f = data(:, 3);

s_ar_goodman = s_a ./ (1 - s_m / s_u);
s_ar_morrow = s_a ./ (1 - s_m / s_fB);

log10_N_f = log10(N_f);
log10_s_ar_goodman = log10(s_ar_goodman);
log10_s_ar_morrow = log10(s_ar_morrow);

disp("2 (a)");
X = [log10_N_f, ones(size(log10_N_f))];
a = (X' * X) \ (X' * log10_s_ar_goodman);
A = 10^(a(2, 1));
B = a(1, 1);
disp("A : ");
disp(A);
disp("B = b : ");
disp(B);

figure(1);
x = 3.5:0.1:6.5;
y = a(1, 1) * x + a(2, 1);
plot(log10_N_f(1:6, 1), log10_s_ar_goodman(1:6, 1), 'kx', log10_N_f(7:13, 1), log10_s_ar_goodman(7:13, 1), 'ko', x, y, 'k');
title("2 (a) Goodman");
xlabel("log_1_0 N_f");
ylabel("log_1_0 \sigma_a_r");

disp("2 (b)");
X = [log10_N_f, ones(size(log10_N_f))];
a = (X' * X) \ (X' * log10_s_ar_morrow);
A = 10^(a(2, 1));
B = a(1, 1);
disp("A : ");
disp(A);
disp("B = b : ");
disp(B);

figure(2);
x = 3.5:0.1:6.5;
y = a(1, 1) * x + a(2, 1);
plot(log10_N_f(1:6, 1), log10_s_ar_morrow(1:6, 1), 'kx', log10_N_f(7:13, 1), log10_s_ar_morrow(7:13, 1), 'ko', x, y, 'k');
title("2 (b) Morrow");
xlabel("log_1_0 N_f");
ylabel("log_1_0 \sigma_a_r");
