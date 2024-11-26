clear
clc
addpath('./solver');
addpath('./quality_assess');
addpath('./tensor_toolbox');

% Load data
datapath = '/MATLAB Drive/NPAlgos/NN''s/NSWTLD/';
addpath(datapath);
matfile = [datapath, 'demodata.mat'];
data = load(matfile);
Omega = data.Omega;

% Load and process image
img_path = [datapath, 'IMG-140935-0001.png'];
hyperspectral_band = im2double(imread(img_path));
E = hyperspectral_band(:,:,1);
E_resized = imresize(E, [200, 200]);
E = repmat(E_resized, [1, 1, size(Omega, 3)]);

% Parameters
Par = struct();
Par.lambda = 0.1;  % Regularization parameter for nuclear norm
Par.Iter = 50;  % Number of iterations

% Call the SNN algorithm 
[Res, Par] = SNN_Completion(E, Omega, Par);  % This is a custom SNN completion function

% Quality assessment
[PSNR, SSIM, SAM] = HSIQA(Res * 255, E * 255);
fprintf('SNN PSNR: %2.3f, SSIM: %2.4f, SAM: %2.3f \n', PSNR, SSIM, SAM);
