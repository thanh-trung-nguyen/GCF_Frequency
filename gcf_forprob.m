function [U] = gcf_forprob(CoefVec,b,WaveNumber,x0)
% function GCF_FORPROB 
% solve the forward problem of the globally convex functional approach for
% the coefficient identification problem of the 1d wave equation in frequency domain. 
% Precisely, it approximates the solution u(x,k) of the following problem
% u''(x,k)  + k^2 c(x) u(x,k) = -delta(x - x0), x in (-infty,infty)
% with the assumption that c(x) = 1 for x outside of (0,b). The solution
% also satisfies the radiation conditions at +- infinity. 
% 
% Method: use the Lippman-Schwinger integral equation: 
% u(x,k) = exp(-ik|x-x0|)/2ik + k^2\int_0^b [exp(-ik|x-y|)/2ik *(c(y)-1) u(y,k) ]dy
% then discretize this integral equation to obtain a linear system. 
%
% INPUT: 
%   CoefVec: a vector of the coefficient c(x). The number of element of
%   this vector is the number of discretization points.
%   b: the upper bound of the interval of x, x in (0,b)
%   WaveNumber: a vector of wavenumbers. 
%   x0: the location of the source, in this setting, we assume that x0 < 0.
% OUTPUT
% U: a matrix of the solution. Each row is the solution at a fixed
% frequency. 
% -------------------------------------------------------------------------
% @Thanh Nguyen, 2016. 



N = length(CoefVec); % number of grid points
Nwn = length(WaveNumber); % number of wave numbers

if size(CoefVec,2) > 1
    CoefVec = CoefVec(:); % convert to column vector if needed
end

h = b/(N-1); % discretization step size
x = linspace(0,b,N); % this is a row vector

U = zeros(Nwn,N); 
CoefMatrix = zeros(N); 

for m = 1:Nwn
    wn = WaveNumber(m); 

    RHS = exp(-1i*wn*abs(x - x0))/2/1i/wn; % right hand side vector
    RHS = RHS(:); % convert to column vector
    
    for n = 1:N
        CoefMatrix(n,:) = exp(-1i*wn*abs(x - x(n))).*(CoefVec.' - 1);
    end
    
    CoefMatrix = eye(N) - wn*h/2/1i*CoefMatrix;
    
    U(m,:) = CoefMatrix\RHS; 
end


