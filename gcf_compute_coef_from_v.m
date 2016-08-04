function Coef = gcf_compute_coef_from_v(v,wn,dx)
% compute the coefficient c(x) from the function v(x) using the equation:
% -c(x) = dv/dx + wn^2*(v(x))^2
% Since v is complex, and c is real, we only take the real part
% 
% Input: v: a vector of values of v(x) at grid point at the given wavenumber
%        which is provided in the second input parameter
%        wn: wavenumber
%        dx: spatial grid size, to approximate the derivative of v(x). 
% 
% Output: Coef: a ROW vector of coefficient values of the same size as v, 
%         note that the last value of Coef is assigned to be 1 (from the
%         theory)
% % =======================================================================



Coef = ones(1,length(v));
Coef(1:end-1) = -real(((v(2:end)-v(1:end-1))/dx + wn^2*v(1:end-1).^2));


