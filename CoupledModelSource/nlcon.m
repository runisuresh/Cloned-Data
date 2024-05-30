function [c,ceq] = nlcon(x)
ceq = [];
c = sin(x(1)) - cos(x(2));
end