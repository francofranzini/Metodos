// metodo_jacobi: Dada la matriz cuadrada A, el vector de TI b, una aproximacion
// inicial x0 y un epsilon para frenar las iteraciones, la funcion devuelve
// el vector solucion x del sistema Ax = b. Esto lo hace mediante el metodo de
// Jacobi.

function x = metodo_jacobi(A, b, x0, eps)

    [filasA, colsA] = size(A);
    [filasB, colsB] = size(b);
    [filasX0, colsX0] = size(x0);

    if filasA <> colsA then // Chequeos de tamaño.
        error("metodo_jacobi: La matriz A debe ser cuadrada.");
        abort;
    end
    if (filasB <> colsA) | (colsB <> 1) then
        error("metodo_jacobi: Tamaño de vector b incorrecto.");
        abort;
    end
    if (filasX0 <> colsA) | (colsX0 <> 1) then
        error("metodo_jacobi: Tamaño de vector x0 incorrecto.");
        abort;
    end

    n = filasA; // Inicializacion.
    x = x0;
    xk = x;
    cont = 0;

    for i = 1:n // Calculo la primer iteracion.
        suma = 0;
        for j = 1:n
            if i<>j then
                suma = suma + A(i,j) * xk(j);
            end
        end
        x(i) = (1/A(i,i))*(b(i)-suma)
    end
    cont = 1;

    while norm(x-xk) > eps // Limite de iteraciones.
        xk = x;
        for i = 1:n
            suma = 0;
            for j = 1:n
                if i<>j then
                    suma = suma + A(i,j) * xk(j);
                end
            end
            x(i) = (1/A(i,i))*(b(i)-suma);
        end
        cont = cont + 1;
    end

    disp(cont);
endfunction
