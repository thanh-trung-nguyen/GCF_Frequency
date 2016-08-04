function f = laguerre_integral_new(n,x,x0)

% calculate the Laguerre function for i=0, 1, 2, 3, ... n. 
% n: the largest index of the Laguerre function to be calculated (start
% from L_0)
% x: a column vector of points where f is calculated
% output: f: a matrix, the first column is L_0(x), the second column is L_1(x), ...
%

x = x(:); % convert to a column vector. 

x = x - x0;
Nx = length(x); 
if (n == 0)
    f = ones(Nx,1);
elseif (n == 1)
    f = [ones(Nx,1), 1-x];
else
%     f = 1/n*((2*n-1-x)*Laguerre_function(n-1,x) - (n-1)*Laguerre_function(n-2,x));
    f = zeros(Nx,n+1);
    
    f(:,1) = 1; 
    f(:,2) = 1 - x;
    for k = 1:n-1;
        f(:,k+2) = ((2*k+1-x).*f(:,k+1) - k*f(:,k))/(k+1); 
    end
end
f = f.*(exp(-x/2)*ones(1,n+1));


% % use the explicit formula: 
% f = laguerre_explicit(n,x,x0);

