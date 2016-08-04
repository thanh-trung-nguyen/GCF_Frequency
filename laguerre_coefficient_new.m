function qn = laguerre_coefficient_new(Nq,q,s,s0)
% % compute the coefficients in the expansion of a function w.r.t. the Laguerre basis in L2: 
% qn(x) = int_{s0}^infty q(s) f_i(s) ds, where f_i(s) = L_i(s-s0), i = 0,
% 1, ..., n
% q can be a COLUMN vector or a matrix whose columns are functions of s.
% s: a vector of frequency, s0: minimum frequency 
% Nq: the largest Laguerre's index
% here we assume that q(s) = O(1/s^2) for s large, and we approximate q(s)
% as follows: 
%   q(s) = C/s^2, where C is calculated as C = q(end)*s(end)^2. 
%



if size(s,1)==1
    s = s(:);  % s must be a column vector
end

smax = 1000; 
ds2 = 0.01; 
ds = s(2) - s(1);

Ns = length(s); 

s2 = (s(end)+ds2:ds2:smax)';

f = laguerre_new(Nq,[s;s2],s0);
f = f'; % for the following computation: 

qn = f(:,1:Ns)*q(:,:)*ds + f(:,Ns+1:end)*((1./s2.^2)*(q(Ns,:)*s(Ns)^2))*ds2;
