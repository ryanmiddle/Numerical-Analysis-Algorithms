%% Simpson's 1/3 Rule Numerical Integration Function
function [I] = Simpson(x,y)
% Ryan Middle's Simpson's 1/3 function 
% This funciton uses the Simpson's 1/3 rule to numerically evaluate the
% integral of output values y with respect to x. If the amount of data
% points is odd, the trapezoidal technique will be used to evaluate the
% final segment of data points.
% Inputs:
% x: the independant variable vector 
% y: the dependant variable vector (with respect to the x vector)
% Outputs:
% I: the numerical integral
if nargin ~= 2
    error('This function requires exactly two inputs')
end
if (size(x) ~= size(y))
    error('The size of each vector must be equal')
end
if length(x) ~= length(y)
    error('Input vectors must be same length')
end
% Determine if inputs are vectors
[m, n]= size(x);
if m ~= 1 && n ~= 1
    error('Inputs must be vectors (1*n or n*1 matrix)')
% For function to work with both row and collumn vectors
elseif m == 1
    iter=n; 
else
    iter=m;
end
% check if vector has even amt of values- requires trap rule
R= rem(iter,2);
if R == 0
    warning('The trapezoidal function will be used for the last segment')
end
% check if vector is evenly spaced
s=diff(x);
if min(s) ~= max(s)
    error('Independant variable is not evenly spaced')
end
I=0;
i=1;
while i<iter-1
    I=I+((x(i+2)-x(i))*(y(i)+4*y(i+1)+y(i+2)))/6;
    
    i=i+2;
end
if R == 0 % Trap Rule
        
       I=I+((x(iter)-x(iter-1))*(y(iter)+y(iter-1)))/2;
        
end

        
