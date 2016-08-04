function [c,qn] = gcf_Laguerre_method(qn1,qn2,x,wn0,Nq,lambda,options,Q_InitGuess)
% solve the inverse problem in frequency domain using the globally convex
% functional algorithm. 
% % Input: 
% qn1: Dirichlet data for q at x = 0, a column vector. Complex function
% qn2: Dirichlet data at x = b, a column vector. Complex function
% x: vector of spatial grid points including two ends 0 and b.
% wn0: the smallest wavenumber. 
% Nq: number of Laguerre's functions
% lambda: Carleman weight parameter
% options: optimization parameters
% Q_0: initial guess (optional)


if size(qn1,2) > 1
    qn1 = qn1.';
end
if size(qn2,2) > 1
    qn2 = qn2.';
end
if size(x,1) > 1
    x = x'; % x is a row vector.
end


h = x(2)-x(1);
Nx = length(x) - 2; % Nx is the number of unknown grid points (excluding the end points with known boundary conditions)


% compute Fkmn and Gkn:
filenameF = sprintf('%s%d%s%3.2f%s','globconvex_coef_F_Nq_',Nq,'_Smin_',wn0,'.mat');
if ~exist(filenameF,'file')
    F = globconvex_coef_F(Nq,wn0);
    eval([' save ' filenameF ' F']);
else
    eval(['load ' filenameF]);
end


% Minimize the objective function:

% for storing the result:
qn = zeros(Nq,Nx+2); % the reconstruction functions qn (complex)
qn(:,1) = qn1;       % the given boundary condition at x = 0
qn(:,Nx+2) = qn2;    % the given boundary condition at x = b


% The initial guess:   
if nargin < 8
    Q_InitGuess = qn2*ones(1,Nx);
end

Q_InitGuess = ComplexMatrix2RealVector(Q_InitGuess); 


CarWF = exp(-lambda*x(1:Nx+1)); % Carleman weight function


% % Check the gradient calculated by the adjoint method:
% dv = 0.01; v = 8:dv:10; Nv = length(v);
% Objfun = 0*v; Grad = zeros(1,Nv); 
% Idx = 1;
% for idk = 1:Nv;
%     InitGuess = [Q_InitGuess(1:Idx-1); v(idk); Q_InitGuess(Idx+1:end)];
%     [Objfun(idk),grad] = objfun_fw(InitGuess,RegPar);
%     Grad(idk) = grad(Idx);    
% end
% figure; plot(v,Objfun); title('Objective function');
% figure; plot(v,Grad(1:Nv)); hold on; plot(v(2:end-1),(Objfun(3:end)-Objfun(1:end-2))/2/dv,'--r');
% legend('Adjoint method','FD approximation'); grid on;



% Solve the optimization problem:
options = optimset(options,'GradObj','on','display','iter','Algorithm','sqp');
%     [Sol,fval] = fmincon(@(Q)objfun_fw(Q,RegPar),Q_InitGuess,[],[],[],[],lb,ub,[],options);
Sol = fminunc(@(Q)objfun_fw(Q),Q_InitGuess,options);

qn(:,2:Nx+1) = RealVector2ComplexMatrix(Sol,Nq,Nx); % take only the last iteration
    
% compute the coefficient c from qn: 
v = compute_v_from_laguerre_coef(qn,wn0,wn0);
c = gcf_compute_coef_from_v(v,wn0,h);

% =========================================================================

    % objective function & gradient:
    function [J,gradJ] = objfun_fw(Q)
   
        Q = [qn1, RealVector2ComplexMatrix(Q,Nq,Nx), qn2];  % add the boundary condition to Q. 
        Scaling = 1; 
       

        J = 0; % objective function
        JR = zeros(Nq,Nx+1); JI = JR;  % store these two matrices for calculating the gradient
        
        
        for m = 1:Nx+1
            for n = 1:Nq
                [SumR,SumI] = sum_F(F,Q(:,m),n);
                JR(n,m) = (real(Q(n,m+1)) - real(Q(n,m)))/h + SumR; 
                JI(n,m) = (imag(Q(n,m+1)) - imag(Q(n,m)))/h + SumI;
                
                J = J + (JR(n,m)^2 + JI(n,m)^2)*CarWF(m);
            end            
        end
        J = h*J*Scaling;

        % gradient:         
        if nargout > 1
            gradR = zeros(Nq,Nx); % gradient of J w.r.t. the real part of Q
            gradI = gradR; % gradient of J w.r.t. the imaginary part of Q

            for m = 1:Nx+1
                for n = 1:Nq                
                    % gradient: 
                    gradJR_QR = zeros(Nq,Nx); % gradient of J(m,n) w.r.t. the real part of Q
                    gradJR_QI = gradJR_QR; % gradient of J(m,n) w.r.t. the imaginary part of Q
                    gradJI_QI = gradJR_QR; % gradient of J(m,n) w.r.t. the imaginary part of Q
                    gradJI_QR = gradJR_QR; % gradient of J(m,n) w.r.t. the imaginary part of Q
                   
                    if (m > 1)             
                        for j = 1:Nq
                            [SumR,SumI] = sum_F_grad(F,Q(:,m),n,j);

                            gradJR_QR(j,m-1) =  SumR; 
                            gradJR_QI(j,m-1) = -SumI; 
                            gradJI_QI(j,m-1) =  SumR; 
                            gradJI_QR(j,m-1) =  SumI; 
                        end
                        gradJR_QR(n,m-1) = gradJR_QR(n,m-1) - 1/h; 
                        gradJI_QI(n,m-1) = gradJI_QI(n,m-1) - 1/h;
                    end
                    if (m <= Nx)
                        gradJR_QR(n,m) = 1/h;
                        gradJI_QI(n,m) = 1/h;
                    end
                    
                    gradR = gradR + 2*CarWF(m)*JR(n,m)*gradJR_QR + JI(n,m)*gradJI_QR; 
                    gradI = gradI + 2*CarWF(m)*JR(n,m)*gradJR_QI + JI(n,m)*gradJI_QI; 

                end 
            end
            % collapse the matrices of gradient and concatenate to a vector: 
            gradJ = h*[gradR(:); gradI(:)]*Scaling;

        end        
    end

    
    function [SumR,SumI] = sum_F(F,Q,n)
        SumR = 0;
        SumI = 0; 
        for j = 1:Nq
            for l = 1:Nq
                SumR = SumR + F(n,j,l)*(real(Q(j))*real(Q(l)) - imag(Q(j))*imag(Q(l)));
                SumI = SumI + F(n,j,l)*(real(Q(j))*imag(Q(l)) + imag(Q(j))*real(Q(l)));
            end
        end
        
    end
 
   function [SumR,SumI] = sum_F_grad(F,Q,n,j)
        SumR = 0;
        SumI = 0; 
        for l = 1:Nq
            SumR = SumR + (F(n,j,l) + F(n,l,j))*real(Q(l));
            SumI = SumI + (F(n,j,l) + F(n,l,j))*imag(Q(l));
        end
        
    end

    % convert a vector to a complex matrix. The vector 
    function Matrix = RealVector2ComplexMatrix(V,nrow,ncol)
        % V is a real vector containing the real part (first half) and
        % imaginary part (second half). 
        
        N = length(V);
        if mod(N,2)~= 0 
            error('In function Vector2ComplexMatrix: length of input vector must be even'); 
        end
        N = N/2; 
        if nrow*ncol ~= N
            error('In function Vector2ComplexMatrix: number of rows and columns must be consistent with the input vector');
        end
        
        RealV = V(1:N); ImagV = V(N+1:2*N);
        Matrix = reshape(RealV,nrow,ncol) + 1i*reshape(ImagV,nrow,ncol);
    end

    function Vector = ComplexMatrix2RealVector(M)
        % M is a matrix of complex numbers. V is a column vector.
        RealM = real(M); 
        ImagM = imag(M); 
        Vector = [RealM(:); ImagM(:)];

    end

        


end

