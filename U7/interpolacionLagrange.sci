// Interpolación Lagrange
function y = Lk(x,k)
    [Xn,Xm] = size(x)
    r = [x(1:k-1) x(k+1:Xm)]
    p = poly(r,"x","roots")
    pk = horner(p,x(k))
    y = p / pk
endfunction

function z = interpolLagrange(x,y)
    [Xn,Xm] = size(x)
    pol = 0
    for k = 1:Xm
        pol = pol + (Lk(x,k)*y(k))
    end
    z = pol
endfunction


function z = error(grado, cota_phi, cota_f)
  z = (cota_phi/prod(grado+1))*cota_f
endfunction

//LINEAL
x0 = [0 0.2]
y0 = [1.0 1.2214]

//CUBICA
x1 = [0 0.2 0.4 0.6]
y1 = [1.0 1.2214 1.4918 1.8221]

lineal = interpolLagrange(x0, y0)
cubica = interpolLagrange(x1, y1)

rango = [-2:0.01:2]

plot(rango,horner(lineal,rango),"r")
plot(rango,horner(cubica,rango),"b")
plot(rango,exp(rango),"g")
a=gca();
a.x_location = "origin";
a.y_location = "origin";
h1 = legend("Lineal","Cubico","e^x")
/*
  Error en el caso lineal es:
    er_lineal(x) = (x-0)*(x-0.2)/2*exp(c_x)''
 
  exp(x)'' = exp(x), tomo el extremo derecho 0.2
  er_lineal(1/3) =  ((1/3)-0)*((1/3)-0.2)/2*exp(0.2) 
                = 0.0271423
*/

//phi = (x - x0)(x-x1)...(x-xn)


//phi <= (ver graficamente)
phi_lineal = poly(x0', 'x')
phi_cubico = poly(x1', 'x')
plot(rango, phi_lineal, "k")
plot(rango, phi_cubico, "k")
e_lineal = error(1, 0.5,1.24)
e_cubico = error(3, 0.02, 1.83)
disp('Error Lineal')
disp(e_lineal)
disp('Error Cubico')
disp(e_cubico)


disp('Error real lineal y cubico')
disp(abs(horner(lineal, 1/3) - 1.395612425))
disp(abs(horner(cubica, 1/3) - 1.395612425))
