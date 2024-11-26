function [Res_NPWTLD, Par] = NL_Denoising(N_Img, Omega, Par)
    % Define the mean filter size directly
    meanFilterSize = [3, 3]; 
    % original Par.rmean with meanFilterSize in the imfilter function
    Xmean = imfilter(N_Img, fspecial('average', meanFilterSize), 'corr', 'replicate', 'same');

    Res_NPWTLD = N_Img;  
end
