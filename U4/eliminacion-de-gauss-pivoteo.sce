// Metodo de eliminacion de gauss con pivoteo parcial.

function [T, g] = eliminacion_gauss_pivoteo(A, b)
    [filA, colA] = size(A);
    [filb, colb] = size(b);
    if (filA <> colA) then
        error("eliminacion_gauss_pivoteo: Matriz no cuadrada");
        abort;
    end
    if (filA <> filb) | (colb <> 1) then
        error("eliminacion_gauss_pivoteo: Tamanio incorrecto de b");
        abort;
    end

    n = filA;
    a = [A b]; // Matriz aumentada

    for k = 1 : n - 1 // Ecuacion k
        kpivot = k;
        amax = abs(a(k,k));  // Pivoteo
        for i = k + 1 : n // Busco el maximo, itero filas.
            if abs(a(i,k)) > amax then
                kpivot = i;
                amax = abs(a(i,k));
            end;
        end;

        if kpivot <> k then
            temp = a(kpivot,:); // Guardo la fila del pivote
            a(kpivot,:) = a(k,:); // En la fila del pivote pongo la fila k actual
            a(k,:) = temp; // En la fila actual pongo el pivote.
        end

        for i = k + 1: n // Fila i
            mult = a(i, k) / a(k, k);
            a(i,k) = 0;
            for j = k + 1: n + 1 // Columna j, incluyendo al vector b.
                a(i,j) = a(i,j) - mult * a(k,j);
            end 
        end
    end;

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

function x = resuelve_sist_elim_gauss_piv(A, b)
    [T, g] = eliminacion_gauss_pivoteo(A, b);
    x = sustitucion_regresiva_superior(T, g);
endfunction
