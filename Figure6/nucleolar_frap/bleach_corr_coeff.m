function bleach_slope = bleach_corr_coeff(post_stack)
%%bleach_corr removed background from stack and corrects for photobleaching
%
%   input :
%       post_stack : A 3D matrixing containing the post-bleach timelapse.
%       Third dimension is time.
%
%   ouput :
%       bleach_slope : A scalar indicating how many intensity units are
%       lost on average per aquisition per pixel.


%% Use mulithresh for mean background subtraction
fluor_means = zeros([size(post_stack,3),1]);
bg_means = zeros([size(post_stack,3),1]);
for t = 1:size(post_stack, 3)
    frame = post_stack(:,:,t);
    bg_threshs = multithresh(frame, 2);
    bg_means(t) = mean(frame(frame < bg_threshs(1)));
    fg_mean = mean(frame(frame > bg_threshs(2)));
    fluor_means(t) = fg_mean - bg_means(t);
end
coeffs = polyfit(1:size(post_stack,3), fluor_means,1);
bleach_slope = coeffs(1);