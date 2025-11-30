// Ejercicio 1 de la Práctica 8
clear()

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

// - Ejercicio i) y ii) - //

disp("a) f(x) = ln(x) y [a,b] = [1,2].")

int_trap = Trapecio(1,2,log)
int_simp = Simpson(1,2,log)
int_real = intg(1,2,log)

disp("Utilizando la Regla del Trapecio se tiene: "+string(int_trap))
disp("Utilizando la Regla de Simpson se tiene: "+string(int_simp))
disp("El valor real (calculado con la función intg) es: "+string(int_real))

// Ejercicio: Calcular el error absoluto y/o relativo.

disp("Error absoluto Trapecio: "+string(int_real - int_trap))
disp("Error absoluto Simpson: "+string(int_real - int_simp))
disp("Error relativo Trapecio: " + string((int_real-int_trap)/2-1))

disp("b) f(x) = x^1/3(x) y [a,b] = [0,0.1].")

function y = f(x)
    y = x^(1/3)
endfunction

int_trap = Trapecio(0,0.1,f)
int_simp = Simpson(0,0.1,f)
int_real = intg(0,0.1,f)

disp("Utilizando la Regla del Trapecio se tiene: "+string(int_trap))
disp("Utilizando la Regla de Simpson se tiene: "+string(int_simp))
disp("El valor real (calculado con la función intg) es: "+string(int_real))

// Ejercicio: Calcular el error absoluto y/o relativo.

disp("Error absoluto Trapecio: "+string(int_real - int_trap))
disp("Error absoluto Simpson: "+string(int_real - int_simp))
disp("Error relativo Trapecio: " + string((int_real-int_trap)/0.1))
disp("Error relativo Simpson: " + string((int_real-int_simp)/0.1))


