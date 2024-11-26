% main
clear
clc
addpath('./solver');
addpath('./quality_assess');
addpath('./tensor_toolbox');
datapath = 'D:\×ÀÃæ\1';
addpath(datapath);

dataname = 'demodata';
load([dataname,'.mat']);

% parameter -- nonlocal similarity
Par.win =   20;   % Non-local patch searching window
Par.delta     =   0;  % Parameter between each iter
Par.Constant         =  2 * sqrt(2);   % Constant num for the weight vector
Par.Innerloop =   2;   % InnerLoop Num of between re-blockmatching
Par.ps       =   20;   % Patch size  best(changed ,the last is 10) 
Par.step        =   5; %  
Par.Iter          =   1; % total iter numbers
Par.nlsp  =  70; % number of K nearest patches

% parameter -- WST
Par.astep = 0.5; % tune
Par.P = [4,1,2,3]; % tune
Par.R1 = [2,2,2,2]; % for Isize [R1,R2,R3,R4]
Par.R2 = [9,9,2,20]; % for Isize [R1,R2,R3,R4] [6,6,2,6] 0.6726 [6,6,6,3] 0.5733
% [6,6,2,3] 0.5821 
% [6,3,2,6] 0.6324
% [6,6,2,10] 0.7199
% [6,6,2,15] 0.7441
% [9,9,2,20] 0.8034 21.719
% [6,9,2,6]  0.6777
% [9,9,2,10] 0.7460
% parameter -- weigth determination
Par.rmean = 5; % window of median filter
Par.gamma = 5; % gamma for partioning Gaussian noise and impulse noise according to noise level
Par.cimpulse = 0.1; % weigth of impulse noise 0.001
Par.eta = 0.05; % for ASTHOSVD
Par.walpha = 0.15; % for linear combination of global and local information

% Initialization
Par.I = data; % clean images
Par.nim = T; % noisy images
clear T data

t1 = tic;
[Res, Par] = NL_Denoising(Par.nim, Omega, Par);
time = toc(t1);
Res(Res>1) = 1;
Res(Res<0) = 0;
% [psnrvec, ssimvec] = Assessment(Res, Par.I);
[PSNR, SSIM, SAM] = HSIQA(Res*255, Par.I*255);
fprintf('PSNR: %2.3f, SSIM: %2.4f, SAM: %2.3f, time: %2.2f \n',PSNR, SSIM, SAM, time);


