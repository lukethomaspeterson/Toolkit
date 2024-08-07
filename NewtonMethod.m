function [x,f,div] = NewtonMethod(fs,dfs,x0,tol,maxIter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NewtonMethod uses the Newton-Raphson method to attempt to find a single %
% root of a function iteratively given some initial guess.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:  fs      = symbolic function                                    %
%          dfs     = symbolic function, derivative of fs                  %
%          x0      = initial guess, dim. must match inputs of fs          %
%          tol     = tolerance for convergence (optional input)           %
%          maxIter = maximum # of iterations before quitting              %
% Outputs: x       = location of the root of f                            %
%          f       = fs(x), function evaluated at the root                %
%          div     = 1 if alogirthm diverged, 0 otherwise                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

div = 0; 

% Optional Inputs Check 
if (nargin < 4 || isempty(tol))
    tol = 1e-15; % set tolerance if none specified 
end

if (nargin < 5 || isempty(maxIter))
    maxIter = 200; % set maximum # of iteration if none are specified 
end

% Newton-Raphson Method
x = x0; % set initial condition 
for i = 1:maxIter 
    
    % Newton-Raphson iteration 
    x = x - fs(x)/dfs(x); 
    f = fs(x); 
    
    % Check Error 
    if abs(f) < tol
        break
    end
end

% Check if the algorithm diverged 
if abs(f) > abs(fs(x0))
    x = x0; 
    f = fs(x0);
    div = 1; % return div = 1
    warning('Algorithm diverged. Returning initial guess.')
end 

end