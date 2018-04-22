%% False Position Algorithm
function falsePosition(func,xl,xu,es,maxiter)
% this function uses the false position method to estimate the zeros of a 
%given funtion. The variables func, xl, xu, es, and iter are required in 
%order for this function to run. if es and iter are not inputed,they will 
%default to 0.0001 error and 200 iterations.
% variable names and meaning:
% inputs:
% func: the function whose roots are trying to be solved for
% xl: independent variable on the lower side of the root 
% xu: independent variable on the upper side of the root
% es: desired relative error (Defaults to 0.0001%)
% maxiter: number of iterations desired (Defaults to 200)
% outputs:
% root: estimated root location
% fx: the function evaluated at the root location
% ea: approximate relative error (%)
% iter: how many iterations were performed
% note that the user must already have his function set with symbolic
% variables

% variable check
if nargin <3 || nargin>5
    error('This function needs at least three inputs and no more than five')
elseif nargin <4 % min amt of variables allowed
    es=0.0001; % default relative error as a percent
    maxiter=200; % default iterations
elseif nargin <5 % user defined relative error but not iterations
    maxiter=200; % default iteratons
else % user defined every variable
end

if func(xl) == 0 % lower limit was the root
    disp('Congrats your lower limit guess is the root!')
elseif func(xu) == 0 % upper limit was the root
    disp('Congrats your upper limit guess is the root!')
else % The initial guesses are not the root of the function
end
iter=0; % no iterations yet so iter is 0
ea=100; % approx relative error is 100% no iters occured yet
root=xl; % initial guess for root is xl, in order to define approx relative error
if func(xl)*func(xu)>0 % did not bind the zero of a function
    error('Your x values chosen are not bounds to a root of this function')
else
while iter<maxiter && es<ea % continue to approximate until error and iteration criteria are met
lastroot=root; % last root used to define the approx relative error
root=xu-(func(xu)*(xl-xu))/(func(xl)-func(xu)); % false position equaiton
ea=abs((root-lastroot)/(root))*100; % calculation of approx relative error as a percent
if func(root)*func(xl)<0 % checks if next bounds need to be between lower and root guess
    xu=root;
else
    xl=root; % the next bound is between the upper and root guess
end
iter=iter+1; % increases current iteration by 1
end

end
fx=func(root);
% display all outputs
display(ea)
display(fx)
display(root)
display(iter)
end
