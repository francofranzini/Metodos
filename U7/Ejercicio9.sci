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

function [p,err] = MinCuad_pol(A,b)
    // Entrada: b = vectores 1xn
    // Salida: p = polinomio de mínimos cuadrados; err = vector de errores (eps = Ax-b)
     [w,a] = gausselimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     err = A*w-b'
endfunction

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

function y = f(x)
    y = 1./(1 + x.^2)
endfunction


for n = [2, 4, 6, 10, 14]
    x = [-5:10/n:5]
    y = f(x)
    
    A = A_mc(x, n)
    [p, e] = MinCuad_pol(A, y)
    disp('-----------------------------')
    disp(p)
    disp(sum(abs(e.^2)))
    disp('-----------------------------')
    intervalo = [-5:0.01:5]
    p_intervalo = horner(p, intervalo)
    f_intervalo = f(intervalo)
    dif_intervalo = f_intervalo - p_intervalo
    dif_intervalo_acotado = dif_intervalo
    //d(intervalo) = f(intervalo) - horner(p, intervalo)
    //plot([-5:0.01:5], dif_intervalo)
    plot(intervalo, dif_intervalo, 'p')
end



