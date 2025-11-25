// Factorizacion LU con pivoteo

function [P, L, U] = fact_LU_pivoteo(A)
    [filasA, colsA] = size(A);
    if filasA <> colsA then
        error("fact_LU_pivoteo: La matriz A debe ser cuadrada.");
        abort;
    end

    n = filasA;
    U = A;
    L = eye(n, n);
    P = eye(n, n);

    for k = 1 : n - 1
        indicePivote = k;
        pivMax = abs(U(k,k));
        for i = k + 1 : n
            pivAnalizar = abs(U(i,k));
            if pivAnalizar > pivMax then
                pivMax = pivAnalizar;
                indicePivote = i;
            end
        end

        if indicePivote <> k then
            tempU = U(indicePivote,:); // Guardo la fila del pivote
            U(indicePivote,:) = U(k,:); // En la fila del pivote pongo la fila k actual
            U(k,:) = tempU; // En la fila actual pongo el pivote.

            if k > 1 then
                tempL = L(indicePivote,1: k-1); // En la matriz L cambio las
                L(indicePivote,1:k-1) = L(k,1:k-1); // filas debajo de la diagonal.
                L(k,1:k-1) = tempL;
            end
            
            tempP = P(indicePivote,:);
            P(indicePivote,:) = P(k,:);
            P(k,:) = tempP;
        end
        
        for j = k + 1: n
            L(j, k) = U(j, k) / U(k, k);
            U(j, k : n) = U(j, k : n) - L(j, k) * U(k, k : n);
            U(j, k) = 0;
        end
    end
endfunction

function x = sustitucion_regresiva_superior(A, b)
    [filA, colA] = size(A);
    [filb, colb] = size(b);
    if (filA <> colA) then
        error("sustitucion_regresiva_superior: Matriz no cuadrada");
        abort;
    end
    if (filA <> filb) | (colb <> 1) then
        error("sustitucion_regresiva_superior: Tamanio incorrecto de b");
        abort;
    end

    n = filA;
    for i = n : -1 : 1
        suma = 0;
        if i <> n then
            for j = i+1:n
                suma = suma + A(i,j) * x(j);
            end
        end
        x(i) = (1/A(i,i))*(b(i) - suma);
    end
endfunction


function x = sustitucion_progresiva_inferior(A, b)
    [filA, colA] = size(A);
    [filb, colb] = size(b);
    if (filA <> colA) then
        error("sustitucion_regresiva_inferior: Matriz no cuadrada");
        abort;
    end
    if (filA <> filb) | (colb <> 1) then
        error("sustitucion_regresiva_inferior: Tamanio incorrecto de b");
        abort;
    end

    n = filA;
    for i = 1 : n
        suma = 0;
        if i <> 1 then
            for j = 1:i-1
                suma = suma + A(i,j) * x(j);
            end
        end
        x(i) = (1/A(i,i))*(b(i) - suma);
    end
endfunction

A = [1 1 1; 1 1 0; 1 0 1];
D = diag(diag(A));
I = eye(3,3);
L = zeros(3, 3);
U = zeros(3,3);
for k = 1:3
    L(k,1:k) = A(k,1:k);    
    U(1:k, k) = A(1:k,k);
end
L = L - D;
U = U - D;

N = D;

T = I - N^-1*A;


disp(max(abs(spec(T)))); 

N = I - (L+D)^-1*A;

T = I - N^-1*A;


disp(max(abs(spec(T))));    

