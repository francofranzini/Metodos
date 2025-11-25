// Metodo de eliminacion de Gauss.

function [T, g] = eliminacion_de_gauss(A,b)
    [filA, colA] = size(A);
    [filb, colb] = size(b);
    if (filA <> colA) then
        error("eliminacion_de_gauss: Matriz no cuadrada");
        abort;
    end
    if (filA <> filb) | (colb <> 1) then
        error("eliminacion_de_gauss: Tamanio incorrecto de b");
        abort;
    end

    n = filA;
    a = [A b];
    for k = 1: n - 1 // Ecuacion k
        for i = k + 1: n // Fila i
            mult = a(i, k) / a(k, k);
            a(i,k) = 0;
            for j = k + 1: n + 1 // Columna j, incluyendo al vector b.
                a(i,j) = a(i,j) - mult * a(k,j);
            end 
        end
    end

    T = a(1 : n, 1 : n)
    g = a(: , n + 1);
endfunction

//-----------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------

function x = resuelve_sist_elim_gauss(A, b)
    [T, g] = eliminacion_de_gauss(A, b);
    x = sustitucion_regresiva_superior(T, g);
endfunction
