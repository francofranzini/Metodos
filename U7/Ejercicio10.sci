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

function w = roots_Cheby(n)
    // Entrada: n = grad del polinomio de Chebyshev
    // Salida: w = vector con las raices del polinomio de Chebyshev
    for i=0:(n-1) do
        w(i+1)=cos((2*i+1)*%pi/(2*n))
    end
endfunction

function w = Cheby(x,n)
    // Entrada: n = número natural; x = número real
    // Salida: Polinomio de Chebyshev de grado n evaluado en x
    if n==0 then
        w = 1
    elseif n==1 then
        w = x
    elseif n==2 then
        w = 2*x.^2-1
    else
        w = 2*x.*Cheby(x,n-1)-Cheby(x,n-2)
    end
endfunction

x = [-1:1]

//Calculamos el polinomio de Chebyshev en el intervalo
//[-1, 1] de grado 4


//No necesitamos el polinomio, sino saber sus raices
//y salen con la siguiente funcion

r = roots_Cheby(4)
//Ahora, para obtener un polinomio que interpole a la funcion f cometiendo el minimo error posible, lo hacemos con los valores de las raices del polinomio de chebyshev, por que? 
/*
    Nosotros buscamos el polinomio que minimice el maximo error cometido por el polinomio, pero obtenerlo significa un problema de optimizacion re costoso. Recordemos que la formula del error esta dada por 
f(x) - pn(x) = (phi_n(x)/(n+1)!)*f^n+1 (Psi(x))
Donde phi_n es el polinomio de raices de nodos de interpolacion, polinomio de grado n, que minimiza su maximo, minimizando la formula del error lo maximo posible.
*/
disp(r)
A = A_mc(r', 2)
[p, e] = MinCuad_pol(A, exp(r)')
disp(p)
e_interval = exp([-1:0.01:1])
p_interval = horner(p, [-1:0.01:1])
err_interval = e_interval - p_interval

disp(sum(abs(e)))
plot([-1:0.01:1],err_interval, "k")
