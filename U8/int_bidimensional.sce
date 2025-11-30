// Ej 5
function v = TnExt(f,a,b,c,d)
    h = (b-a)*(d-c) / 4
    v = h * (f(c,a)+f(c,b)+f(d,a)+f(d,b))
endfunction

function f = f1(x,y)
    f = sin(x+y)
endfunction

// falta dividirlo en intervalos flaco (Skinny)
funcprot(0)

function v = Tn(fx,a,b,n)
    h = (b-a)/n;
    v = (fx(a)/2)+ (fx(b)/2);
    for i = 1 : n-1
        xi = a+i*h;
        v = v + fx(xi);
    end
    v = h*v;
endfunction

function v = DobleTn(f,a,b,c,d,n,m) // nab, mcd, (c, d funciones)
    h = (b-a)/n;
    deff("z=fxa(y)","z=f("+string(a)+",y)")
    deff("z=fxb(y)","z=f("+string(b)+",y)")
    v = (Tn(fxa,c(a),d(a),m) + Tn(fxb,c(b),d(b),m))/2
    
    for i = 1:n-1
        xi = a + i*h
        deff("z=fxi(y)","z = f(" +string(xi)+ ",y)")
        v = v + Tn(fxi,c(xi),d(xi),m)
    end
    v = h * v
endfunction

function v=DobleTn_opt(f,a,b,c,d,n,m)
    h = (b-a)/n
    
    function t = Tn_fijo_x(x)
      function fy = f_fijo(y)
        fy = f(x,y);
      endfunction
      t = Tn(f_fijo,c(x),d(x),m);    
    endfunction
      
    
    v = (Tn_fijo_x(a) + Tn_fijo_x(b)) / 2;
    for i = 1: n-1
        xi = a+i*h;
        v = v + Tn_fijo_x(xi);
    end
    v = h*v;
endfunction


function v=Sn(fxx,x0,x1,n)
    h = (x1-x0)/n
    v = (fxx(x0)+fxx(x1))
    for j = 1:1:n-1
        xi = x0+j*h
        if modulo(j,2) == 1 then
            v = v + 4*fxx(xi)
        else
            v = v + 2*fxx(xi)
        end
    end
    v = v*h/3
endfunction

function v = DobleSn(f,a,b,c,d,n,m) // nab, mcd, (c, d funciones)
    h = (b-a)/n;
    deff("z=fxa(y)","z=f("+string(a)+",y)")
    deff("z=fxb(y)","z=f("+string(b)+",y)")
    
    v = Sn(fxa,c(a),d(a),m) + Sn(fxb,c(b),d(b),m)
    for i = 1:n-1
        xi = a + i*h
        deff("z=fxi(y)","z = f(" +string(xi)+ ",y)")
        if modulo(i,2) == 1 then
            v = v + 4*Sn(fxi,c(xi),d(xi),m)
        else
            v = v + 2*Sn(fxi,c(xi),d(xi),m)
        end
    end
    v = (h/3) * v
endfunction

function v = DobleSn_opt(f,a,b,c,d,n,m) // nab, mcd, (c, d funciones)
    h = (b-a)/n;
    
    function t = fijar_x(x)
      function fy = fijar_y(y)
        fy = f(x,y);
      endfunction
      t = Sn(fijar_y,c(x),d(x),m);
    endfunction
    
    v = fijar_x(a) + fijar_x(b);
    
    for i = 1:n-1
        xi = a + i*h
        if modulo(i,2) == 1 then
            v = v + 4 * fijar_x(xi);
        else
            v = v + 2 * fijar_x(xi);
        end
    end
    v = (h/3) * v
endfunction


function y = cx(x)
    y = -sqrt(2*x-x^2);
endfunction

function y = dx(x)
    y = sqrt(2*x-x^2);
endfunction

function v = uno(x,y)
    v = 1;
endfunction

//res = DobleTn(uno,0,2,cx,dx,100,100)

function y = c1(x)
    y = 2
endfunction

function y = d1(x)
    y = 3
endfunction

function f = f2(x,y)
    f = %e^((atan(abs(x*y))/y^2)*(%e^y)*sqrt(abs((%e^x)*x)))
endfunction

function f = f3(x,y)
    f = %e^(sin(tan(x*y)))
endfunction

function h = h_func(x, y)
    deff("I1 = integral_interna(t)", "I1 = (t^2 - 5*t) / (t^4 - 2*t^2 - 1)");
    deff("I2 = integral_externa(x)", "I2 = exp(integrate(integral_interna, t, "+string(x)+", 2))");
    h = exp(integrate(integral_externa, 'j', 0, y)) + %pi;
endfunction

function y=zero(x)
         y = 0
endfunction

function z=uno(x,y)
         z = 1
endfunction

function c = cena(x,y)
  c = sin(x+y);
endfunction

// DobleTn(f,a,b,c,d,n,m)
tic(); int1 = DobleSn(cena,0,1,zero,uno,100,100); t1 = toc();
disp(int1);
disp(t1);
tic(); int2 = DobleSn_opt(cena,0,1,zero,uno,100,100); t2 = toc()
disp(int2);
disp(t2);
