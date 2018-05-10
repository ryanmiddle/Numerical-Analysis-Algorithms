function [t,y] = Heun(dydt,tspan,y0,h,es,maxit)
% This algorithm uses Heun's method to numerically solve ordinary
% differential equations. By inputting a differential equation with
% variables y and t, with a domain tspan, as well as an initial condition
% and a specified height. es and maxit are optional conditions and will be
% filled if not specified. These allow Heun's method to iterate until
% iteration or error requirements are met. The output is a data table with
% values y with respect to t, which are then displayed on a graph.
% Inputs:
% dydt- differential equation handle with variables y and t
% tspan- the domain of the function where Heun's method will be applied.
% Note that if tspan and h values do not agree, meaning that the domain is
% not evenly divisible by h, the final approximation will use a step size h
% that is smaller than the one desired by the user.
% y0- initial condition for diffeq at initial t condition given in variable
% tspan
% h- the height between t values 
% es- min approx. relative error for Heun's predictor-corrector. Defaults
% to 0.001
% maxit- maximum amount of iterations allowed for Heun's predictor
% corrector. Defaults to 50
% Outputs:
% t- a vector with quantities for the independant variable, specified by
% tspan and h
% y- a vector with quantities for the dependant variable with respect to t
tic
if nargin <4 || nargin > 6
    error('This algorithim requires a minimum of 4 arguments and a maximum of 6 arguments')
elseif nargin == 4
    es=0.001; % Default relative error
    maxit=50; % Default max iterations
else
    maxit=50;
end
ti=tspan(1); % Min t value/starting condition
tf=tspan(2); % Max t value
t=ti:h:tf; % creates t vector with every point
n=length(t); 
if t(n) ~= tf % checks if final point agrees with step size, adjusts final iter step size
    t(n+1)=tf;
    n=n+1;
end
y=zeros(1,n);
slopes=zeros(1,n);
for i=1:n-1 % Euler's Method to predict new y value
 
    y(i)=y0; 
    
    h=t(i+1)-t(i); % determines step size for each point
    slopes(i)=dydt(t(i),y(i)); % left slope for current iteration is slopes(i)
    y(i+1)=slopes(i)*h+y(i); % prediction of rightmost y value
    iter=1; % initalizes current iteration
    ea=1; % initalized current error
    while (iter<maxit) && (es<abs(ea)) % Heun's Predictor Corrector for y(i+1)
        slopes(i+1)=dydt(t(i+1),y(i+1)); % Righthand slope 
         ynew=y(i)+(h/2)*(slopes(i)+slopes(i+1)); % predictor-corrector
        ea=(ynew-y(i+1))/ynew; % error of current estimate
        iter=iter+1; % iterations done on predictor-corrector
        y(i+1)=ynew;
    end
   y0=ynew; % stores predictor-corrected y value into y0 so it can be placed in the y vector during next iter
 
end
plot(t,y,'m*--')
xlim([ti, tf])
xlabel('t')
ylabel('y')
legend('Heuns method for differential equation dydt')    
  
   
 toc  
end 
