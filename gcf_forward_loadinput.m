function [x,NoiseLevel,WaveNumber,SourceLocation,Coef] = gcf_forward_loadinput(paramfile)
% load the input parameters for the 1D wave equation

% -----------------load the parameters: 
fid = fopen(paramfile);
if fid == -1
    error('The parameter file not found');
end

% spatial grid: 
text = fscanf(fid,'%s',1); 
Xmin = str2double(fscanf(fid,'%s',1));
Xmax = str2double(fscanf(fid,'%s',1)) ;
Nx = round(str2double(fscanf(fid,'%s',1)));

x = linspace(Xmin,Xmax,Nx); 

% location of point source:
text = fscanf(fid,'%s',1);
SourceLocation = str2double(fscanf(fid,'%s',1));


% wavenumber: 
text = fscanf(fid,'%s',1); 
WnMin = str2double(fscanf(fid,'%s',1));
WnMax = str2double(fscanf(fid,'%s',1));
Nfreq = round(str2double(fscanf(fid,'%s',1)));

WaveNumber = linspace(WnMin,WnMax,Nfreq)';


% Noise level:
text = fscanf(fid,'%s',1);
NoiseLevel= str2double(fscanf(fid,'%s',1));

% Number of objects:
text = fscanf(fid,'%s',1); 
NrObjects = round(str2double(fscanf(fid,'%s',1)));

% load the object's coefficients: 
Xmin_obj = zeros(NrObjects,1); 
Xmax_obj = Xmin_obj;
Coef_obj = Xmin_obj;
obj_types = cell(NrObjects,1); % types of objects

for i = 1:NrObjects
    text = fscanf(fid,'%s',1);  %#ok<*NASGU>
    Xmin_obj(i) = str2double(fscanf(fid,'%s',1));
    Xmax_obj(i) = str2double(fscanf(fid,'%s',1));
    Coef_obj(i) = str2double(fscanf(fid,'%s',1));
    obj_types{i} = fscanf(fid,'%s',1); 
end

fclose(fid);


% create the coefficient: 
Coef = ones(1,Nx);

for n = 1:NrObjects
    idx = find(x >= Xmin_obj(n) & x <= Xmax_obj(n));
    if strcmpi(obj_types{i},'gaussian')
        Coef(idx) = Gaussian_coefficient(x(idx),Xmin_obj(1),Coef_obj(1),Xmax_obj(1));
    else
        Coef(idx) = Coef_obj(n);
    end
end

   

