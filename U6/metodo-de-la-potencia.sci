// Metodo de la potencia

function [v, lambda] = metodo_potencia(A, z0, eps, maxIter)
    [filA, colA] = size(A);
    [filZ0, colZ0] = size(z0);
    if (filA <> colA) then
        error("metodo_potencia: Matriz no cuadrada");
        abort;
    end
    if (filA <> filZ0) | (colZ0 <> 1) then
        error("metodo_potencia: Tamanio incorrecto de z0");
        abort;
    end
    n = filA;
    
    i = 1;
    // Calculamos la primera iteracion
    wn = A*z0;
    zn = wn / norm(wn, %inf); 
    
    while (i <= maxIter) & (norm(z0-zn) > eps)
        k = 1;
        while(wn(k) == 0) & (k <= n)
            k = k + 1; // Suponemos que existe siempre una componente no nula.
        end
        if (k == n + 1) then
            error("metodo_potencia: wn se convirtio en el vector nulo.");
            abort;
        end
        wn = A*zn;
        z0 = zn; // Guardo el anterior zn.
        zn = wn / norm(wn, %inf);
        v = zn;
        lambda = wn(k) / z0(k);
        i = i + 1;
    end
    printf("Iteraciones %d\n", i);
endfunction
