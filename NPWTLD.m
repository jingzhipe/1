clear
clc
addpath('./solver');
addpath('./quality_assess');
addpath('./tensor_toolbox');

% Load the dataset
datapath = '/MATLAB Drive/NPAlgos/NN''s/NSWTLD/';
addpath(datapath);
matfile = [datapath, 'demodata.mat'];
data = load(matfile);
Omega = data.Omega;

% Load and process the image
img_path = [datapath, 'IMG-140935-0001.png'];
hyperspectral_band = im2double(imread(img_path));

% Reduce image size for faster testing
E = hyperspectral_band(:,:,1);  % First band (grayscale)
E_resized = imresize(E, [100, 100]);  % Resize for faster computation
E = repmat(E_resized, [1, 1, size(Omega, 3)]);  % Replicate for all bands

% Display sizes of E and Omega to identify mismatches
fprintf('Size of E: %s\n', mat2str(size(E)));
fprintf('Size of Omega: %s\n', mat2str(size(Omega)));

% Check if sizes are equal
if ~isequal(size(E), size(Omega))
    % If sizes don't match, reshape or resize E to match Omega's dimensions
    E = imresize(E, [size(Omega, 1), size(Omega, 2)]);
    fprintf('Reshaped E to: %s\n', mat2str(size(E)));
end

% Re-check sizes after reshaping
fprintf('Size of E after reshaping: %s\n', mat2str(size(E)));
fprintf('Size of Omega: %s\n', mat2str(size(Omega)));

% Initialize Par structure with parameters for NPWTLD
Par.win = 20;                
Par.delta = 0;              
Par.Constant = 2 * sqrt(2);
Par.Innerloop = 1;           
Par.ps = 10;                 
Par.step = 3;                
Par.Iter = 1;                
Par.nlsp = 10;               
Par.astep = 0.5;             
Par.P = [4, 1, 2, 3];        
Par.R1 = [2, 2, 2, 2];       
Par.R2 = [9, 9, 2, 20];      
Par.rmean = 5;               
Par.gamma = 5;               
Par.cimpulse = 0.1;          
Par.eta = 0.05;              
Par.walpha = 0.15;

% Preallocate result matrix Res to avoid dynamic memory allocation
Res = zeros(size(E));

% Display sizes of Res to verify it matches E and Omega
fprintf('Size of Res: %s\n', mat2str(size(Res)));

% Perform element-wise multiplication (ensure sizes match)
for i = 1:size(E, 3)
    Res(:,:,i) = E(:,:,i) .* Omega(:,:,i);
end

% Now run the denoising process
t1 = tic;  % Start timing

% Simplify and reduce complexity inside NL_Denoising as needed
[Res, Par] = NL_Denoising(Res, Omega, Par);

% Measure elapsed time
time = toc(t1);

% Clip results to range [0, 1]
Res(Res > 1) = 1;
Res(Res < 0) = 0;

% Quality assessment (PSNR, SSIM)
[PSNR_NPWTLD, SSIM_NPWTLD, SAM_NPWTLD] = HSIQA(Res * 255, E * 255);
fprintf('NPWTLD PSNR: %2.3f, SSIM: %2.4f, SAM: %2.3f, time: %2.2f \n', PSNR_NPWTLD, SSIM_NPWTLD, SAM_NPWTLD, time);

