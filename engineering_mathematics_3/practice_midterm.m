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

A0 = 100 * (rand(200, 200) - 0.5);
%A0 = eye(10); 
%A0 = [1, 1; 1, 0];
%A0 = [3, 1, 2; -2, -5, -5; 4, 7, 1];
%A0 = [5, -2, 1, 4; -2, 3, 0, 1; 1, 0, 2, -3; 4, 1, -3, 2];

lambda = eigenval(A0);
lambda

%% 프로베니우스 행렬

function A = frobenius(coeffs)
    n = size(coeffs, 2);
    A = zeros(n);
    
    A(1, :) = -coeffs;
    A(2:n, 1:n-1) = eye(n-1);
end

%% 하우스홀더 변환

function A = householder(A0)

    n = size(A0, 1);
    A = A0;
    for i=1:(n-1)
        gamma = sign(A(i+1, i)) * sqrt((A((i+1):n, i)' * A((i+1):n, i)));
        u = zeros(n, 1);

        u((i+1):n) = A((i+1):n, i);
        u(i+1) = u(i+1) + gamma;

        alpha = 1 / (gamma^2 + gamma * A(i+1, i));
        if abs(gamma) < 10^-16
            P = eye(n);
        else
            P = eye(n) - alpha * (u * u');
        end
            

        A = P * A * P;

    end
end

%% 이중 QR 반복법

function lambda = double_qr(A0, lambda0)
    max_iter = 10000;

    nn = size(A0, 2);
    if nn == 0
        lambda = lambda0;
        return;
    end
    if nn == 1
        lambda = [lambda0; A0(1, 1)];
        return;
    end
    if nn == 2
        tr = trace(A0);
        de = det(A0);
        lambda = [lambda0; tr / 2 + sqrt(tr^2 / 4 - de)];
        lambda = [lambda; tr / 2 - sqrt(tr^2 / 4 - de)];
        return;
    end
    if nn > 2
        for iter = 1:max_iter
            if iter == max_iter
                disp("did not converge");
                lambda = lambda0;
                return;
            end
            flag = 0;
            for i = 1:nn-1
                if abs(A0(i+1, i)) < 10^(-14)

                    lambda = double_qr(A0(1:i, 1:i), lambda0);
                    lambda = double_qr(A0(i+1:nn, i+1:nn), lambda);
                    flag = 1;

                    break;
                end
            end
            if flag == 1
                break;
            end
            sk0 = A0(nn-1, nn-1);
            sk1 = A0(nn, nn);

            [Q, R] = qr((A0 - sk0 * eye(nn)) * (A0 - sk1 * eye(nn)));

            A0 = Q' * A0 * Q;
        end
    end
end

%% eigenvalues

function lambda = eigenval(A0)
    lambda = [];
    A_householder = householder(A0);
    lambda = double_qr(A_householder, lambda);
end