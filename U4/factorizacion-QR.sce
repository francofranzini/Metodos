// Factorizacion QR

function [Q, R] = factorizacion_QR(A)
    [m, n] = size(A);
    Q = zeros(m,n);
    R = eye(n,n);

    for k = 1:n
        suma = 0
        for i = 1: k - 1
            suma = suma + (A(:,k)' * Q(:,i)) * Q(:,i);
        end
        vk = norm(A(:,k) - suma);
        Q(:,k) = (A(:,k) - suma) / vk;
        R(k,k) = vk;
        for j = k + 1 :n
            R(k, j) = A(:,j)' * Q(:, k);
        end
    end
endfunction
