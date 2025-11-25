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

x = [0, 1, 2, 2.5]
y = [1, 3, 3, 3]

P_0123 = interpolLagrange(x,y)

disp(horner(P_0123, 3))
