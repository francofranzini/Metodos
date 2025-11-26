// Ejercicio 7 de la Práctica 7
clc()
clear()

// Aproximación polinomial de mínimos cuadrados polinomial para matrices con rango completo
function [p,err] = MinCuad_pol(A,b)
    // Entrada: b = vectores 1xn
    // Salida: p = polinomio de mínimos cuadrados; err = vector de errores (eps = Ax-b)
     [w,a] = gausselimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     err = A*w-b'
endfunction

// Matriz del método de mínimo cuadrados polinomial
function A = A_mc(x,n)
    // Entrada: x,y = vectores 1xn; n = grado de la aproximación
    // Salida: A = matriz del método de mínimo cuadrados
    m = length(x)
    A = ones(m,1)
    //[Phi_1(x1:xn), Phi_2(x1:xn), ..., Phi_n(x1:xn)]
    for j=2:(n+1) do
        A = [A,(x').^(j-1)]
    end
endfunction

// Método de Eliminación Gaussiana con pivoteo parcial
function [x,a] = gausselimPP(A,b)
[nA,mA] = size(A) 
[nb,mb] = size(b)
a = [A b]; // Matriz aumentada
n = nA;    // Tamaño de la matriz
// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    kpivot = k; amax = abs(a(k,k));  //pivoteo
    for i=k+1:n
        if abs(a(i,k))>amax then
            kpivot = i; amax = a(i,k);
        end;
    end;
    temp = a(kpivot,:); a(kpivot,:) = a(k,:);
    a(k,:) = temp
    
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k)
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end
end
// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n)
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k)
    end
    x(i) = (a(i,n+1)-sumk)/a(i,i)
end
endfunction

x = [0, 0.15, 0.31, 0.5, 0.6, 0.75]
y = [1, 1.004, 1.031, 1.117, 1.223, 1.422]

Amc_1 = A_mc(x, 1)
[p1,e1] = MinCuad_pol(Amc_1, y)
//disp(p1)
disp('ERROR 1')
disp(sum(e1.^2))

Amc_2 = A_mc(x, 2)
[p2, e2] = MinCuad_pol(Amc_2, y)
disp('ERROR 2')
disp(sum(e2.^2))

Amc_3 = A_mc(x, 3)
[p3, e3] = MinCuad_pol(Amc_3, y)

disp('ERROR 3')
disp(sum(e3.^2))



