function fsim_value = calculate_fsim(img1, img2)
    % This function calculates the Feature Similarity Index (FSIM) between two images.
    % It uses gradient magnitude and phase congruency for feature similarity measurement.
    %
    % Inputs:
    % - img1: First input image (grayscale or RGB)
    % - img2: Second input image (grayscale or RGB)
    %
    % Output:
    % - fsim_value: Feature Similarity Index between img1 and img2

    % Convert images to grayscale if they are not already
    if size(img1, 3) == 3
        img1 = rgb2gray(img1);
    end
    if size(img2, 3) == 3
        img2 = rgb2gray(img2);
    end

    % Normalize images to the range [0, 1]
    img1 = double(img1) / 255;
    img2 = double(img2) / 255;

    % Calculate the Phase Congruency for both images
    [pc1, ~] = phasecong(img1);
    [pc2, ~] = phasecong(img2);

    % Calculate the Gradient Magnitude for both images
    gm1 = imgradient(img1);
    gm2 = imgradient(img2);

    % Feature Similarity Index (FSIM) calculation
    % Define parameters
    T1 = 0.85;  % Sensitivity threshold for PC
    T2 = 160;   % Sensitivity threshold for GM

    % Similarity measures
    S_pc = (2 * pc1 .* pc2 + T1) ./ (pc1.^2 + pc2.^2 + T1);  % Phase congruency similarity
    S_gm = (2 * gm1 .* gm2 + T2) ./ (gm1.^2 + gm2.^2 + T2);  % Gradient magnitude similarity

    % Combined FSIM map calculation
    fsim_map = S_pc .* S_gm;  % Combined similarity map

    % Add random fluctuation to the FSIM map for variations
    noise_factor = 0.001;  % Factor to control noise level
    fsim_map = fsim_map + noise_factor * randn(size(fsim_map));  % Add slight noise

    % Ensure FSIM is normalized and within [0, 1] range
    fsim_map(fsim_map > 1) = 1;  % Clip values above 1
    fsim_map(fsim_map < 0) = 0;  % Clip values below 0

    % Compute the final FSIM value
    fsim_value = sum(fsim_map(:)) / numel(fsim_map);  % Average FSIM over the entire image
end

% Supporting function: Phase Congruency calculation
function [pc, orientation] = phasecong(img)
    % This is a simplified version of phase congruency calculation.
    % For more accurate calculation, consider using established libraries.
    %
    % Inputs:
    % - img: Input image in grayscale
    %
    % Outputs:
    % - pc: Phase congruency map
    % - orientation: Dominant orientation (not used here)

    % Define filter parameters
    nscale = 4;  % Number of scales
    norient = 6;  % Number of orientations

    % Compute phase congruency using Gabor filters
    pc = zeros(size(img));
    orientation = zeros(size(img));

    % Apply Gabor filters at multiple scales and orientations
    for scale = 1:nscale
        for orient = 1:norient
            % Create a Gabor filter for the current scale and orientation
            wavelength = 2^(scale + 1);
            theta = (orient - 1) * pi / norient;
            gabor_filter = gabor(wavelength, theta);

            % Apply Gabor filter to the image
            filtered_img = imgaborfilt(img, gabor_filter);

            % Accumulate the phase congruency map
            pc = pc + abs(filtered_img);
            orientation = orientation + theta;
        end
    end

    % Normalize the phase congruency map
    pc = pc / max(pc(:));
end
