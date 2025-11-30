function v = Tn(f, x ,a,b,n)
    //x = valor fijo en eje x para integrar en y
    h = (b-a)/n;
    v = (f(x, a)/2)+ (f(x, b)/2);
    for i = 1 : n-1
        xi = a+i*h;
        v = v + f(x, xi);
    end
    v = h*v;
endfunction


function v=DobleTn(f,a,b,c,d,n,m)
    h = (b-a)/n
    v = Tn(f,a,c(a),d(a),m)/2+Tn(f,b,c(b),d(b),m)/2
    for i=1:n-1
        xi = a+i*h;
        v = v + Tn(f, xi,c(xi),d(xi),m);
    end
    v = h*v;
endfunction

function y=zero(x)
         y= 0
endfunction

function y=uno(x)
         y = 1
endfunction

function y = dos(x)
    y = 2
endfunction

function z = f(x,y)
    z = sin(x+y)
endfunction

dobletn = DobleTn(f,0, 2, zero, uno, 2,2)

disp(dobletn)
