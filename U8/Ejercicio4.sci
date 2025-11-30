funcprot(0)
// Regla del Trapecio
function w = Trapecio(a,b,f)
    // Entrada: a,b = extremos de integración; f = función de Scilab
    // Salida: w = aproximación de la integral de f en [a,b] por la Regla del Trapecio
    w = ((b-a)/2)*(f(a)+f(b))
endfunction

// Regla de Simpson
function w = Simpson(a,b,f)
    // Entrada: a,b = extremos de integración; f = función de Scilab
    // Salida: w = aproximación de la integral de f en [a,b] por la Regla de Simpson
    h = (b-a)/2
    w = (h/3)*(f(a)+4*f(a+h)+f(b))
endfunction

function w = compuestoTrapecio(a, b, f, n)
    h = (b-a)/n
    w = (1/2)*f(a) + (1/2)*f(b)
    for j = 1:n-1 do
        w = w + f(a+j*h)
    end
    w = h*w
    
endfunction

function w = compuestoSimpson(a,b,f,n)
    h = (b-a)/n
    w = 0
    if n == 2 then
        w = Simpson(a, b, f)
    else
        //w = f(x0) + f(xn)
        w = f(a) + f(b)
        for j = 1:n-1
            if modulo(j, 2) == 0 then
               w = w + 2*f(a+j*h) 
            else
               w = w + 4*f(a+j*h)
            end
        end
        //w = f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + ... + f(xn)
        w = (h/3)*w
    end
endfunction

function y = errorCompuestoTrap(a, b, n,cota_f)
    h = (b-a)/n
    y = -((h^2)*(b-a)/12) *cota_f
endfunction

function y = errorCompuestoSimp(a, b, n, cota_f)
    h = (b-a)/n
    y = -((h^4)*(b-a)/180)*cota_f
endfunction

function y = f(x)
    y = 1/(x+1)
endfunction
a = 0
b = 1.5
n = 10
I = 0.916290737
int_trap = compuestoTrapecio(a,b,f,n)
int_simp = compuestoSimpson(a,b,f,n)

disp("Utilizando la Regla del Trapecio se tiene: "+string(int_trap))
disp("Utilizando la Regla de Simpson se tiene: "+string(int_simp))
disp("El valor real (calculado con la función intg) es: "+string(I))

disp('El error estimado obtenido por el compuesto trapecio es: ' + string(abs(errorCompuestoTrap(a,b,n,2))))
disp('El error estimado obtenido por el compuesto simpson es: ' + string(abs(errorCompuestoSimp(a,b,n,24))))

disp('El error absoluto obtenido por el compuesto trapecio es: ' + string(abs(int_trap - I)))
disp('El error absoluto obtenido por el compuesto simpson es: ' + string(abs(int_simp-I)))
