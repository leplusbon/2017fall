%% °ÅµìÁ¦°ö¹ý

max_iter = 100;

A = [
    3, 9, 6;
    -2, -5, -5;
    4, 7, 17
    ];

x = [1; 1; 1];
lambda = 0;

for i = 1:max_iter
    
    y = A * x;
    lambda = max(y);
    
    x = y / lambda;
    
    disp(x');
    disp(lambda);
    disp("");
    
end

%% ¿ª°ÅµìÁ¦°ö¹ý

max_iter = 100;

A = [
    3, 1, 2;
    -2, -5, -5;
    4, 7, 1
    ];

x = [1; 1; 1];
lambda = 0;

for i = 1:max_iter
    
    y = A \ x;
    
    lambda = max(y);
    
    x = y / lambda;
    
    disp(x');
    disp(1/lambda);
    disp("");
    
end


%% ÇÏ¿ì½ºÈ¦´õ º¯È¯

A0 = 50 * (rand(100, 100) - 0.5);

n = size(A0, 1);

A = A0;

for i=1:(n-1)
    gamma = sign(A(i+1, i)) * sqrt((A((i+1):n, i)' * A((i+1):n, i)));
    u = zeros(n, 1);
    
    u((i+1):n) = A((i+1):n, i);
    u(i+1) = u(i+1) + gamma;
    
    alpha = 1 / (gamma^2 + gamma * A(i+1, i));
    
    P = eye(n) - alpha * (u * u');
    
    A = P * A * P;
    
end

for index = 1:n
    
    nn = 
    for i=1:max_iter
        A_small = A((n-1):n, (n-1):n);
        sk0 = A(n-1, n-1);
        sk1 = A(n, n);

        [Q, R] = qr((A - sk0 * eye(n)) * (A - sk1 * eye(n)));

        A = Q' * A * Q;
        
        if i == max_iter
            disp("did not converge");
            break;
        end
        if abs(A(n, n-1)) < 10^(-14)
            break;
        end

    end
    if i == max_iter
        break;
    end
    lambda(index) = A(n, n);
    disp(lambda(index);
    
    A = A(1:n-1, 1:n-1);
    
end