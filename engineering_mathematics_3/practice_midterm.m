%% 거듭제곱법

max_iter = 100;

A = [
    3, 9, 6;
    -2, -5, -5;
    4, 7, 17
    ];

x = [1; 1; 1];
lambda = 0;

for i = 1:max_iter
    
    lambda_old = lambda;
    y = A * x;
    lambda = max(y);
    
    if abs(lambda_old - lambda) < 10^-14
        break;
    end
    
    x = y / lambda;
    
    fprintf("iterations : %d\n\n", i);
    fprintf("eigenvalue : ");
    disp(lambda);
    fprintf("eigenvector : ");
    disp(x');
    
end

%% 역거듭제곱법

max_iter = 100;

A = [
    3, 1, 2;
    -2, -5, -5;
    4, 7, 1
    ];

x = [1; 1; 1];
lambda = 0;

for i = 1:max_iter
    
    lambda_old = lambda;
    y = A \ x;
    lambda = max(y);
    
    if abs(lambda_old - lambda) < 10^-14
        break;
    end
    
    x = y / lambda;
    
    fprintf("iterations : %d\n\n", i);
    fprintf("eigenvalue : ");
    disp(1 / lambda);
    fprintf("eigenvector : ");
    disp(x');
    
end


%% 하우스홀더 변환, 이중 QR 분해

clear;

max_iter = 5000;

%A0 = 100 * (rand(100, 100) - 0.5);
%A0 = [3, 1, 2; -2, -5, -5; 4, 7, 1];

A0 = frobenius([-8000.002, 16]);

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

lambda = zeros(n, 1);

AA = A;
index = 1;
for j = 1:size(A, 2)
    
    if index > size(A, 2)
        break;
    end
    
    nn = n - index + 1;
    if nn == 1
        lambda(index) = AA(nn, nn);
        disp(lambda(index));
        index = index + 1;
        AA = [];
        small = 0;

        break;
    end
    for i=1:max_iter

        sk0 = AA(nn-1, nn-1);
        sk1 = AA(nn, nn);

        [Q, R] = qr((AA - sk0 * eye(nn)) * (AA - sk1 * eye(nn)));

        AA = Q' * AA * Q;
        
        if i == max_iter
            disp("did not converge");
            break;
        end
        if nn > 1
            if abs(AA(nn, nn-1)) < 10^(-14)
                small = 1;
                break;
            end
        end
        if nn > 2
            if abs(AA(nn-1, nn-2)) < 10^(-14)
                small = 2;
                break;
            end
        end
        if nn == 2
            small = 2;
            break;
        end

    end
    if i == max_iter
        break;
    end
    if small == 1
        lambda(index) = AA(nn, nn);
        AA = AA(1:nn-1, 1:nn-1);
        disp(lambda(index));
        index = index + 1;
    end
    if small == 2
        A_small = AA(nn-1:nn, nn-1:nn);
        tr = trace(A_small);
        de = det(A_small);
        lambda(index) = tr / 2 + sqrt(tr^2 / 4 - de);
        lambda(index + 1) = tr / 2 - sqrt(tr^2 / 4 - de);
        AA = AA(1:nn-2, 1:nn-2);
        disp(lambda(index));
        disp(lambda(index + 1));
        index = index + 2;
    end
    
end

%% 프로베니우스 행렬

function A = frobenius(coeffs)
    n = size(coeffs, 2);
    A = zeros(n);
    
    A(1, :) = -coeffs;
    A(2:n, 1:n-1) = eye(n-1);
end

