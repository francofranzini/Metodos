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
  z = (cota_phi/prod([1:grado+1]))*cota_f
  disp(z)
endfunction

x = [-1, 1, 2, 4]

y = [2, 4, -1, 37]

P3 = interpolLagrange(x, y)

disp(horner(P3,0))

phi = poly(x', 'x')
disp(phi)

disp(max(abs(horner(phi, [-1, 0.01,4]))))

plot([-1,0.01,4], phi)

error_P3 = error(3, abs(horner(phi, 0)),33.6)

disp(error_P3)
