function [L, U] = dootlittle(A)
    [filasA, colsA] = size(A);
    if filasA <> colsA then
        error("dootlittle: La matriz A debe ser cuadrada.");
        abort;
    end

    n = filasA;
    U = zeros(n, n);
    L = eye(n, n);

    for i = 1:n
        for j = 1:n
            if i <= j then
                sumatoria = 0;
                for k = 1: i-1
                    sumatoria = sumatoria + L(i,k)* U(k,j);
                end
                U(i,j) = A(i,j) - sumatoria;
            else
                sumatoria = 0;
                for k = 1: j-1
                    sumatoria = sumatoria + L(i,k)* U(k,j);
                end
                L(i,j) = (A(i,j) - sumatoria) / U(j,j);
            end
        end
    end
endfunction
