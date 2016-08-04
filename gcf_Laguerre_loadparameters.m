function [FolderName,X,WaveNumber,numLagFunc,CalermanWeight,...
            CoefTruncThreshold,CoefLowerBound,CoefUpperBound,MaxIter]...
                                   = gcf_Laguerre_loadparameters(paramfile)
% load the input parameters for the Laguerre's methods for the CIP:

% -----------------load the parameters: 
fid = fopen(paramfile);
if fid == -1
    error('The parameter file not found');
end

% folder name: 
text = fscanf(fid,'%s',1);
FolderName = fscanf(fid,'%s',1);

% spatial domain: 
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
Xmin = str2double(fscanf(fid,'%s',1));
Xmax = str2double(fscanf(fid,'%s',1)) ;
Nx = round(str2double(fscanf(fid,'%s',1)));
X = linspace(Xmin,Xmax,Nx);


% interval of wavenumbers:
text = fscanf(fid,'%s',1); 
WnMin = str2double(fscanf(fid,'%s',1));
WnMax = str2double(fscanf(fid,'%s',1));
Nfreq = round(str2double(fscanf(fid,'%s',1)));

WaveNumber = linspace(WnMin,WnMax,Nfreq)';


% number of Laguerre's functions: 
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
numLagFunc = round(str2double(fscanf(fid,'%s',1)));

% Coefficients of the Carleman weight function:
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
CalermanWeight = str2double(fscanf(fid,'%s',1));

% Coefficient truncation threshold:
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
CoefTruncThreshold = str2double(fscanf(fid,'%s',1));

% lowerbound and upper bound
text = fscanf(fid,'%s',1);  %#ok<*NASGU>
CoefLowerBound = str2double(fscanf(fid,'%s',1));
CoefUpperBound = str2double(fscanf(fid,'%s',1));


% maximum number of iterations
text = fscanf(fid,'%s',1); 
MaxIter = round(str2double(fscanf(fid,'%s',1)));
    


fclose(fid);

