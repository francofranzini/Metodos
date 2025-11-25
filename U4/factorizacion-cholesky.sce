// Factorizacion de Cholesky

function U = cholesky(A)
    [filA, colA] = size(A);
    if filA <> colA then
        error("cholesky: La matriz A debe ser cuadrada");
        abort;
    end
    
    n = filA;
    eps = 10e-8; // Epsilon para saber si A no es definida positiva.
    U = zeros(n,n);

    for k = 1 : n
        // Calculo de U(k,k)
        if k == 1 then
            t = A(1,1);
        end
        if k <> 1 then
            sumatoria = 0;
            for i = 1 : k - 1
                sumatoria = sumatoria + U(i,k)**2;
            end
            t = A(k,k) - sumatoria;
        end
        
        if t < eps then
            error("cholesky: Matriz A no definida positiva");
            abort;
        end
        U(k,k) = sqrt(t);
        
        // Calculo del resto de la fila/columna.
        for j = k + 1:n
            if k == 1 then
                U(1,j) = A(1,j) / U(1,1);
            end
            if k <> 1 then
                sumatoria = 0;
                for i = 1 : k - 1
                    sumatoria = sumatoria + U(i,k)*U(i,j);
                end
                U(k,j) = (A(k,j) - sumatoria) / U(k,k);
            end
        end
    end
endfunction
