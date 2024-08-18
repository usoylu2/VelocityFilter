function [peak_locations, I, J] = get_centroids_coef(frameGray)

peakmap1 = normxcorr2_same(frameGray);
peakmap1(peakmap1 <20000) = 0;
figure
imagesc(peakmap1)
title('Crosscorrelation in Passband');
ylim([125 205]);
%figure;
%imagesc(abs(hilbert(peakmap1)))
%title('Crosscorrelation after Enveleope detection');
%peakmap1 = abs(hilbert(peakmap1));
%peakmap1(peakmap1 < 3*max(max(peakmap1))/4) = 0; 

%peakmap1(peakmap1 <16000) = 0;
%peakmap1(peakmap1 < 0.12) = 0; %for linear case
% % [x, y] = find(coef_map == max(max(coef_map)))
peak_locations = imregionalmax(peakmap1,8);    % local max     % size(BW) - size(frameGray);
% BW = imregionalmax(frameGray,8);
peakvalues = [];
counter =1;
for i =1:numel(peak_locations)
    if (peak_locations(i) ==1)
        peakvalues(counter) = peakmap1(i);
        counter = counter +1;
    end
end
peak_locations = peak_locations .* peakmap1;
peak_locations(peak_locations <prctile(peakvalues,80)) = 0;
f = find(peak_locations);
[I, J] = ind2sub(size(peak_locations),f);

function coef_map = normxcorr2_same(frame)      % normalized cross correlation, output the same size....
load('PSFrf.mat')
dz = 3.08e-05;
dx = 3.08e-05;
Nx = 1238;
Nz = 975;
psf = data(end - ceil(9.5e-3/ dz)+2:end- ceil(7.5e-3/ dz)+1, Nx/2 - ceil(1e-3/dx)+1 : Nx/2 +ceil(1e-3/dx)-1);
size(psf);
psf = flip(psf,1);
psf = flip(psf,2);
coef_map= conv2(frame,psf);
%coef_map = normxcorr2(psf, frame);
coef_map = coef_map(((size(psf, 1) - 1 )/ 2) + 1 : size(coef_map,1) - ( (size(psf, 1)-1)/ 2), ...
    ((size(psf, 2)-1)/ 2) + 1 : size(coef_map,2) - ((size(psf, 2)-1)/ 2));
% clip the edge.
