function F = practice_fun(x)
k = 1:10;
F = 2 + 2*k - exp(k*x(1)) - 2*exp(2*k*(x(2)^2));
end