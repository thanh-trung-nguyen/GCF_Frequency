function test_gcf_freq_sim(parameterfile_forprob, parameterfile_Laguerre,StepToRun)
% test the globally convex functional algorithm in frequency domain.
% input: 
%   parameterfile_forprob: parameter file for forward problem. This is
%   used ONLY for simulated data
%   parameterfile_Laguerre: parameter file for the GCF algorithm. 
%   StepToRun: step to run, choose between 1 and 2, 1 if you want to simulate the data, 2 if
%   you only run the inverse algorithm. If StepToRun >= 2, the first
%   parameter file is ignored, so it can be provided as an empty variable
%   []
% OUTPUT: several output data files and figures in the folder specified in
% the parameter file for the GCF algorithm. 
% =======================================================================


% parameterfile_forprob = 'parameter_forprob.dat';
% parameterfile_Laguerre = 'parameter_inversion_Laguerre.dat';
% parameterfile_Tikhonov = 'parameter_inversion_Tikhonov.dat';


if nargin < 3
    error('Not enough input parameters. Must be 3');
end

if StepToRun == 1
    %=== Load parameters for forward solver: for simulated data only
    [x,NoiseLevel,WaveNumberFor,SourceLocation,CoefExact] ...
                            = gcf_forward_loadinput(parameterfile_forprob);

    Nx = length(x); % number of spatial grid points, including the endpoints
    dx = x(2)-x(1); % spatial grid step size

end

if StepToRun <= 2
    % ==== Load the input parameter for the GCF method:
 
    [FolderName,X_Lag,WaveNumber,NrLagFunc,CalermanWeight,CoefTruncThreshold,...
        CoefLowerBound,CoefUpperBound,MaxIterLag] = gcf_Laguerre_loadparameters(parameterfile_Laguerre);

    NxLag = length(X_Lag); % number of grid points for the coefficient. 
    dxLag = X_Lag(2)-X_Lag(1);
    
    Nfreq = length(WaveNumber); % number of frequencies 
    ds = WaveNumber(2)-WaveNumber(1); % step size in the wave number (~frequency)

    % check the consistency of the wavenumbers in the two parameter files:
    if exist('WaveNumberFor','var') && max(abs(WaveNumber - WaveNumberFor)) > 3*eps
        error('The wavenumbers in the two parameters files must be the same');
    end    
        
    % create the directory to store files: 
    if ~exist(FolderName,'dir')  % create the folder if not exist.
        eval(['mkdir ' FolderName]);
    end   

end

% ==== Step 1: Simulate the data: 
if StepToRun == 1
    disp('Solving the forward problem and simulating the measured data.')
    disp('If you want to run the inversion for this data set next time, choose StepToRun > 1');
    
    % solve the forward problem: 
    exactSol = gcf_forprob(CoefExact,x(Nx),WaveNumber,SourceLocation);
    
    % create the data for inverse problem: 
    DirichletData = exactSol(:,1); 
    
    % calculate the Neumann data:
    C = DirichletData  - exp(-1i*WaveNumber*abs(SourceLocation))/2/1i./WaveNumber; 
    C = C.*(1 + NoiseLevel*abs(C)); % add random multiplicative noise to the scattered wave
    
    DirichletData = exp(-1i*WaveNumber*abs(SourceLocation))/2/1i./WaveNumber + C; 
    NeumannData = -exp(-1i*WaveNumber*abs(SourceLocation))/2 + 1i*WaveNumber.*C;
    
    % save the Dirichlet and Neumann data. % the first column is the real
    % part, the second column is the imaginary part. 
    dlmwrite([FolderName,'/DirichletData.dat'],[real(DirichletData), imag(DirichletData)],'delimiter',' ');
    dlmwrite([FolderName,'/NeumannData.dat'],[real(NeumannData), imag(NeumannData)],'delimiter',' ');

    % --- the following part is only for checking the calculation: 
    figure(1); set(gca,'fontsize',15);
    hold off; plot(WaveNumber,real(NeumannData));
    hold on; plot(WaveNumber,real((exactSol(:,2) - exactSol(:,1))/dx),'--r'); hold off;
    title('Compare the calculation of Neumann data'); 
    legend('using closed formula','FDM approximation');
    xlabel('x'); ylabel('real part of the Neumann data'); 
    
    % calculating v, q and qn:   
    v = zeros(Nfreq,Nx);     % function v = u_x/k^2/u
    v(:,1:end-1) = (exactSol(:,2:end) - exactSol(:,1:end-1))/dx./exactSol(:,1:end-1)./(WaveNumber.^2*ones(1,Nx-1)); 
    v(:,Nx) = v(:,Nx-1); % v(x) = -1i/k = constant for x=a at which coefficient = 1 for x >= a.
    
    q = (v(2:Nfreq,:) - v(1:Nfreq-1,:))/ds;
    qn = laguerre_coefficient_new(NrLagFunc-1,q,WaveNumber(1:end-1),WaveNumber(1)); % "exact" Laguerre's coefficients
    
    % ---- compute the coefficient epsilon(x) from the Laguerre's coefficient: 
    qnLag = zeros(NrLagFunc,length(X_Lag));
    for k = 1:NrLagFunc
        qnLag(k,:) = linearinterpolation(qn(k,:),x,X_Lag);
    end
    
    qLag = compute_q_from_laguerre_coef(qnLag,WaveNumber(1:end-1),WaveNumber(1));
    vLag = compute_v_from_laguerre_coef(qnLag,WaveNumber(1:end-1),WaveNumber(1));
    CoefLagExact = gcf_compute_coef_from_v(vLag(1,:),WaveNumber(1),dxLag);
    CoefLagExact(CoefLagExact <CoefLowerBound) = CoefLowerBound; CoefLagExact(CoefLagExact > CoefUpperBound) = CoefUpperBound;

    set(gca,'fontsize',15); 
    plot(X_Lag,CoefLagExact,'-b','linewidth',2); 
    hold on; plot(x,CoefExact,'--r','linewidth',2); 
    title('Compare computed and exact coefficients'); 
    xlabel('x'); ylabel('c(x)'); 
    legend('Computed','Exact'); 
    
    
end % end of Step 1.


% === Step 2: Running the Laguerre's method:
if StepToRun <= 2    
    disp('Running the Laguerre method: ');
    
    % ----load the data: 
    DirichletData = dlmread([FolderName,'/DirichletData.dat']);  % Dirichlet data u. This is the total wave
    DirichletData = DirichletData(:,1) + 1i*DirichletData(:,2); 
    NeumannData = dlmread([FolderName,'/NeumannData.dat']); % Neumann data du/dx
    NeumannData = NeumannData(:,1) + 1i*NeumannData(:,2); 

    % compute v, q, and the Laguerre coefficients qn: 
    v0 = NeumannData./WaveNumber.^2./DirichletData; % boundary condition for function v at x = 0.
    q0 = (v0(2:Nfreq) - v0(1:Nfreq-1))/ds;
    qb = 1i./WaveNumber.^2;    
    
     % --- Laguerre's coefficients of the boundary data: 
    qn1 = laguerre_coefficient(0:NrLagFunc-1,q0,WaveNumber(1:end-1),WaveNumber(1)); % boundary data at x = 0
    qn2 = laguerre_coefficient(0:NrLagFunc-1,qb,WaveNumber,WaveNumber(1)); % boundary data at x = b
   
    
    % --- Main part of step 2: Minimizing the Laguerre's functional: 
    options = optimset('MaxIter',MaxIterLag,'TolFun',1e-10);
    [CoefLag,qnr] = gcf_Laguerre_method(qn1,qn2,X_Lag,WaveNumber(1),NrLagFunc,CalermanWeight,options,qnLag(:,2:end-1));
    CoefLag(CoefLag <CoefLowerBound) = CoefLowerBound; CoefLag(CoefLag > CoefUpperBound) = CoefUpperBound;

    plot(X_Lag,CoefLag); title('reconstructed coefficient');
    
    
%     % --- Compute the coefficient in the whole interval, this will be used as the initial guess for the local method
%     Index1 = find(X_FDM > X_Lag(2),1,'first');
%     Index2 = find(X_FDM < X_Lag(end-1),1,'last');
%     CoefLagFull = 0*CoefExact+1;
%     CoefLagFull(Index1:Index2) = linearinterpolation((CoefLag-1)*1.8+1,X_Lag(2:end),X_FDM(Index1:Index2));  
% 
    % --- display the function q and Laguerre's coefficients: 
    figure; set(gca,'fontsize',16);
    for n = 1:NrLagFunc
        plot(X_Lag, qnr(n,:),'linewidth',2); hold on; plot(X_Lag, qnLag(n,:),'--r','linewidth',2); hold off;
        legend('Comuted Laguerre coefficient','Exact Laguerre coefficient');
        xlabel('X'); ylabel('Laguerre coefficient');
%         print('-depsc2',[FolderName,'/Lag_coef_',num2str(n),'.eps']);
    end
% 
%     % display the coefficient of Laguerre's method together with the true
%     % one: TEST ------------------
%     CoefLagTrunc = CoefLagFull; CoefLagTrunc(CoefLagFull < 0.75*max(CoefLagFull)) = 1; 
%     figure; set(gca,'fontsize',16);
%     idx = length(X_FDM):-1:1;
%     hold off; plot(-X_FDM(idx),CoefExact(idx),'-r','linewidth',2); 
%     hold on; plot(-X_FDM(idx),CoefLagFull(idx),'--b','linewidth',2); hold off;
%     hold on; plot(-X_FDM(idx),CoefLagTrunc(idx),'-.k','linewidth',2); hold off; 
%     legend('Exact coefficient','Result of Step 1','Result of Step 1 with 75% truncation');
%     axis([-X_Lag(end) -X_FDM(1) min([min(CoefExact), min(CoefLagFull)])-0.5 max([max(CoefExact), max(CoefLagFull)]+2)]);
%     xlabel('X'); ylabel('c(x)'); grid on;             
%     print('-depsc2',[FolderName,'/Lag_method.eps']);
    

end
%  end of Laguerre's method.








% 
% % parameters for steps 3 and 4: 
% [RegPar_local_L2,RegPar_local_H1,NrRun,NrIter,NrRefine,MaxIter] = loadparameters_ls(parameterfile_Tikhonov);
% options = optimset('MaxIter',MaxIter,'TolFun',1e-6,'Algorithm','sqp'); % options for the fmin
%         
% if step <= 3    
%     disp('Running the local method using Laguerre result: !!!');    
%     
% %     % --- The result of the Laguerre's method is used as  the initial guess
% %     SubGrid = X_Lag; % 
% %     InitGuess = [0; (CoefLag - 1); 0]; 
% 
%       % no mesh refinement: 
%       SubGrid = X_FDM(Index1:Index2); 
%       InitGuess = CoefLagFull(Index1:Index2) - 1;
%    
%     
%     for RunIdx = 1:NrRun   
%         for n = 1:NrIter
%             % lower and upper bounds:     
%             lb = 0*InitGuess-1 + CoefLowerBound; ub = lb + CoefUpperBound;
% 
%             % run the local method:     
%             [CoefLocal,CoefValues,fval,gradf] = waveeq1d_inverse_ls(Ui_tt,InitGuess,SubGrid,DirichletData,IdxMea,X_FDM,dt,lb,ub,options,RegPar_local_L2,RegPar_local_H1);
% 
%             % refine the mesh for the next iteration: uniform refinement
%             SubGrid = linearinter(SubGrid,NrRefine);
%             InitGuess = linearinter(CoefValues,NrRefine);
% 
% %             % Test: refine the mesh at the points of large gradient only:
% %             SubGridOld = SubGrid; 
% %             SubGrid = mesh_refinement_1d(SubGridOld,gradf,0.5); 
% %             InitGuess = linearinterpolation(CoefValues,SubGridOld,SubGrid);            
% %             plot(SubGrid,InitGuess,'*'); hold on; plot(X_FDM,CoefExact); hold off;
% %             pause;
% 
%             % --- display the coefficient: 
%             figure; set(gca,'fontsize',16);
%             idx = length(X_FDM):-1:1;
%             hold off; plot(-X_FDM(idx),CoefExact(idx),'-r','linewidth',2); 
%             hold on; plot(-X_FDM(idx),CoefLocal(idx),'--b','linewidth',2); hold off;
%             if RunIdx == 1
%                 hold on; plot(-X_FDM(idx),CoefLagFull(idx),'-.k','linewidth',2); hold off; 
%                 legend('Exact coefficient','Result of Step 2','Result of Step 1');
%                 axis([-X_Lag(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact), min(CoefLagFull)])-0.5 max([max(CoefLocal), max(CoefExact), max(CoefLagFull)]+2)]);
%             else
%                 legend('Exact coefficient','Result of Step 3');
%                 axis([-X_Lag(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact)])-0.5 max([max(CoefLocal), max(CoefExact)]+2)]);
%             end
%             
%             xlabel('X'); ylabel('c(x)'); grid on;             
%             print('-depsc2',[FolderName,'/Run',num2str(RunIdx),'_Iter',num2str(n),'.eps']);
%             
%         end % finish each run        
%         
%          % ----after each run, we shrink the interval:
%         CoefLocal(abs(CoefLocal-1) < Coef_Truncation_Threshold) = 1;  
%         minvalue = min(CoefLocal); 
%         maxvalue = max(CoefLocal);
%         if abs(X_mea - X_Lag(end)) < EPS % right-hand side measurement                
%             % find the right-most interval in which the coefficient is large 
%             if minvalue >= 1 - Coef_Truncation_Threshold % strong target   
%                 disp('test');
%                 LeftIdx = find(CoefLocal == maxvalue,1,'first'); RightIdx = LeftIdx; 
%             else % weak target
%                 LeftIdx = find(CoefLocal == minvalue,1,'first'); RightIdx = LeftIdx; 
%             end                          
%             while LeftIdx > 1 && abs(CoefLocal(LeftIdx) - 1) > Coef_Truncation_Threshold
%                 LeftIdx = LeftIdx - 1;
%             end
% %             LeftIdx = LeftIdx - 3                
%             
%             SubGrid = X_FDM(LeftIdx:IdxMea);
%             InitGuess = CoefLocal(LeftIdx:IdxMea);
%             NrIter = 1;
%             
%         elseif abs(X_mea - X_Lag(1)) < EPS % left hand measurement:
%             RightIdx = find(abs(CoefLocal-1) > Coef_Truncation_Threshold,1,'first'); LeftIdx = RightIdx;
%             while RightIdx < NrX_FDM && abs(CoefLocal(RightIdx)-1) > Coef_Truncation_Threshold
%                 RightIdx = RightIdx + 1;
%             end
%             %RightIdx = RightIdx + 5;       
%             SubGrid = X_FDM(IdxMea:RightIdx);
%             InitGuess = CoefLocal(IdxMea:RightIdx);
%             NrIter = 1;
%         end           
%         
%     end    
% end % end of the local method
% 
% 
% if step <= 4    
%     disp('Running the local method alone: !!!');
%     
%    
% %     % --- The result of the Laguerre's method is used as  the initial guess
% %     SubGrid = X_Lag; % 
% %     InitGuess = [0; 0*(CoefLag - 1); 0]; 
%       SubGrid = X_FDM(Index1:Index2); 
%       InitGuess = 0*CoefLagFull(Index1:Index2);
% 
%     ParBasis = FEM_basis(SubGrid,X_FDM); % the basis for parametrizing the coefficient
%     CoefInitGuess = FEM_basis_expansion(ParBasis,InitGuess);
%     
%     for RunIdx = 1:NrRun   
%         for n = 1:NrIter
%             % lower and upper bounds:     
%             lb = 0*InitGuess-1 + CoefLowerBound; ub = lb + CoefUpperBound;
% 
%             % run the local method:     
%             [CoefLocal,CoefValues,fval,gradf] = waveeq1d_inverse_ls(Ui_tt,InitGuess,SubGrid,DirichletData,IdxMea,X_FDM,dt,lb,ub,options,RegPar_local_L2,RegPar_local_H1);
% 
%             % refine the mesh for the next iteration: uniform refinement
%             SubGrid = linearinter(SubGrid,NrRefine);
%             InitGuess = linearinter(CoefValues,NrRefine);
% 
% %             % Test: refine the mesh at the points of large gradient only:
% %             SubGridOld = SubGrid; 
% %             SubGrid = mesh_refinement_1d(SubGridOld,gradf,0.85); 
% %             InitGuess = linearinterpolation(CoefValues,SubGridOld,SubGrid);            
%             
%             % --- display the coefficient: 
%             % --- display the coefficient: 
%             figure; set(gca,'fontsize',16);
%             idx = length(X_FDM):-1:1;
%             hold off; plot(-X_FDM(idx),CoefExact(idx),'-r','linewidth',2); 
%             hold on; plot(-X_FDM(idx),CoefLocal(idx),'--b','linewidth',2); hold off;
%             if RunIdx == 1
%                 hold on; plot(-X_FDM(idx),CoefInitGuess(idx),'-.k','linewidth',2); hold off; 
%                 legend('Exact coefficient','The local method','Initial guess');
%                 axis([-X_LagOld(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact), min(CoefInitGuess)])-0.5 max([max(CoefLocal), max(CoefExact), max(CoefInitGuess)]+2)]);
%             else
%                 legend('Exact coefficient','The local method');
%                 axis([-X_LagOld(end) -X_FDM(1) min([min(CoefLocal), min(CoefExact)])-0.5 max([max(CoefLocal), max(CoefExact)]+2)]);
%             end            
%             xlabel('X'); ylabel('c(x)'); grid on;   
% 
%             print('-depsc2',[FolderName,'/Local_Run',num2str(RunIdx),'_Iter',num2str(n),'.eps']);
%             
%         end % finish each run        
%         
%          % ----after each run, we shrink the interval:
%         CoefLocal(abs(CoefLocal-1) < Coef_Truncation_Threshold) = 1;  
%         minvalue = min(CoefLocal); 
%         maxvalue = max(CoefLocal);
%         if abs(X_mea - X_Lag(end)) < EPS % right-hand side measurement                
%             % find the right-most interval in which the coefficient is large 
%             if minvalue >= 1 - Coef_Truncation_Threshold
%                 LeftIdx = find(CoefLocal == maxvalue,1,'first'); RightIdx = LeftIdx; 
%             else 
%                 LeftIdx = find(CoefLocal == minvalue,1,'first'); RightIdx = LeftIdx; 
%             end                          
%             while LeftIdx > 1 && abs(CoefLocal(LeftIdx) - 1) > Coef_Truncation_Threshold
%                 LeftIdx = LeftIdx - 1;
%             end
%             %LeftIdx = LeftIdx - 5;                
%             
%             SubGrid = X_FDM(LeftIdx:IdxMea);
%             InitGuess = CoefLocal(LeftIdx:IdxMea);
%             NrIter = 1;
%             
%         elseif abs(X_mea - X_Lag(1)) < EPS % left hand measurement:
%             RightIdx = find(abs(CoefLocal-1) > Coef_Truncation_Threshold,1,'first'); LeftIdx = RightIdx;
%             while RightIdx < NrX_FDM && abs(CoefLocal(RightIdx)-1) > Coef_Truncation_Threshold
%                 RightIdx = RightIdx + 1;
%             end
%             %RightIdx = RightIdx + 5;       
%             SubGrid = X_FDM(IdxMea:RightIdx);
%             InitGuess = CoefLocal(IdxMea:RightIdx);
%             NrIter = 1;
%         end           
%         
%     end    
% end % end of the local method
% 
% 
% 
% close all;

