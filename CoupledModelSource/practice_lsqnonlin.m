rng default

d = linspace (0,3);

y = exp(-1.3*d) + 0.05*randn(size(d));

% Anonymous function
% Takes a value of the experimental decay rate, r
% Returns a vector of differences:
% The model with that decay rate and the data
fun = @(r)exp(-d*r)-y;

% Arbitrarily choose an initial guess
x0 = 4;

% Find the value of the optimal decay rate
x = lsqnonlin(fun,x0);

plot(d,y,'ko',d,exp(-x*d),'b-')
legend('Data','Best fit')
xlabel('t')
ylabel('exp(-tx)')


% Find a centering b and scaling a that best fit the function:
% aexp(−t)exp(−exp(−(t−b)))

t = linspace(-4,4);

% Normal density at each point t
y = 1/sqrt(2*pi)*exp(-t.^2/2);

% Create a function that evaluates the difference:
% The centered and scaled function from the normal y
% x(1) as the scaling, a
% x(2) as the centering, b
fun = @(x)x(1)*exp(-t).*exp(-exp(-(t-x(2)))) - y;


% Find the optimal fit starting from x0 = [1/2,0]
x0 = [1/2,0];
% Scaling a between 1/2 and 3/2
% Centering b between -1 and 3
lb = [1/2,-1];
ub = [3/2,3];

x = lsqnonlin(fun,x0,lb,ub);

plot(t,y,'r-',t,fun(x)+y,'b-')
xlabel('t')
legend('Normal density','Fitted function')


% Sum from k = 1 to 10:
% 2 + 2*k - exp(k*x(1)) - 2*exp(2*k*(x(2)^2))
% Check practice_fun

% Minimize this function subject to the linear constraint:
% x1 ≤ x2 / 2
% Written as x1 - x2 / 2 ≤ 0
% Ax ≤ b
A = [1 -1/2]; % 1 * x1 + (-1/2) * x2
b = 0;

% Impose bounds
% 0 ≤ x1 ≤ 2
% 0 ≤ x2 ≤ 4
lb = [0 0];
ub = [2 4];

% Start the optimization process from the point x0 = [0.3 0.4]
x0 = [0.3 0.4];

% The problem has no linear equality constraints
Aeq = [];
beq = [];

% Run the optimization
x = lsqnonlin(@practice_fun,x0,lb,ub,A,b,Aeq,beq);

% Nonlinear constraint:
% sin(x1) ≤ cos(x2)

% No linear constraints
A = [];
b = [];

x = lsqnonlin(@practice_fun,x0,lb,ub,A,b,Aeq,beq,@nlcon);


% Suppose that you have:
% observation time data xdata
% observed response data ydata
% Find parameters x(1) and x(2) to fit a model of the form:
% ydata=x(1)exp(x(2)xdata)

xdata = ...
 [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = ...
 [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];

% Create a simple exponential decay model
% Computes a vector of differences between predicted and observed values.
fun = @(x)(x(1)*exp(x(2)*xdata)-ydata)./ydata;

% Fit the model using the starting point x0 = [100,-1]
% First, use the default 'trust-region-reflective' algorithm.
x0 = [100,-1];
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
x = lsqnonlin(fun,x0,[],[],options);

% See if there is any difference using the 'levenberg-marquardt' algorithm
% (There's not)
options.Algorithm = 'levenberg-marquardt';
x = lsqnonlin(fun,x0,[],[],options);

% Plot the solution and the data
plot(xdata,ydata,'ko')
hold on
tlist = linspace(xdata(1),xdata(end));
plot(tlist,x(1)*exp(x(2)*tlist),'b-')
xlabel xdata
ylabel ydata
title('Exponential Fit to Data')
legend('Data','Exponential Fit')
hold off

