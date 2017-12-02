sn = [675, 578, 600, 559, 563, 540;
    14000, 55000, 58000, 61000, 165000, 270000]';
log_sn = log10(sn);
X = [log_sn(:, 2), ones(size(log_sn, 1), 1)];
a = (X' * X) \ (X' * log_sn(:, 1));

A = 10^(a(2, 1));
B = a(1, 1);

x = 4:0.1:5.5;
y = a(1, 1) * x + a(2, 1);

disp("3(a)");
figure(1);
plot(log_sn(:, 2), log_sn(:, 1), 'k.');
title("3(a)");
xlabel("log_1_0(N_f)");
ylabel("log_1_0(\sigma_a)");
disp("approximate A : ");
disp(10^3);
disp("approximate B : ");
disp(-0.07);

disp("3(b)");
figure(2);
plot(log_sn(:, 2), log_sn(:, 1), 'k.', x, y, 'r');
title("3(b)");
xlabel("log_1_0(N_f)");
ylabel("log_1_0(\sigma_a)");
disp("refined A : ");
disp(A);
disp("refined B : ");
disp(B);
disp("sigma_f' : ");
disp(A / 2^B);
disp("b : ");
disp(B);
