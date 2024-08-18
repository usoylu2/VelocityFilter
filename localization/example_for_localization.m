%% Single Bubble

% f0 = 5e6; % Transducer center frequency [Hz]
% c = 1540; % Speed of sound [m/s]
% lambda = c/f0;
% dx = lambda/10;  
% dz = lambda/10;
% Nx = 1238;
% Nz = 975;
% i= 50;
%temp = u(: ,: ,i);
%a =size(temp);
%figure;
%imagesc(log(abs(hilbert(temp))));
%imagesc(data);
%title('psf');
%ylim([150 180]);
%set(gcf, 'Position', get(0, 'Screensize'));
% [~, row1, col1] = get_centroids_coef(temp);
% hold on
% plot(col1, row1, 'r.')

% ylabel(' Depth(mm)');
% xlabel(' Lateral(mm)');
%dat = data(end - ceil(13.5e-3/ dz):end- ceil(3.5e-3/ dz), Nx/2 - ceil(5e-3/dx) : Nx/2 +ceil(5e-3/dx) ,:);

%figure;
%imagesc(abs(hilbert(dat(:,:,i))));
%imagesc( dat(:,:,i));
%set(gcf, 'Position', get(0, 'Screensize'));

%hold on
%plot(col1, row1, 'r.')
%title('b(r,t)');
%ylim([150 180]);

%% Crossing Vessels
% % clc

f0 = 5e6; % Transducer center frequency [Hz]
c = 1540; % Speed of sound [m/s]
lambda = c/f0;
dx = lambda/10;  
dz = lambda/10;
Nx = 1238;
Nz = 975;
i= 75;
temp = u(: ,: ,i);
a =size(temp);
% figure;
% imagesc(abs(hilbert(temp)));
% title('\phi(r,t)');
%ylim([150 180]);
%set(gcf, 'Position', get(0, 'Screensize'));
[~, row1, col1] = get_centroids_coef(temp);
%columnmodified = (col1-1) * dx*1e3 -5.0204  ;  
%rowmodified = (row1-1) * dz*1e3 +16.5 ; 
figure;
%imagesc(linspace(0,(a(2)-1)*dx*1e3,a(2)) -(a(2)-1)*dx*1e3/2,linspace(0,(a(1)-1)*dz*1e3 ,a(1)) +16.5,abs(hilbert(temp)));
imagesc(abs(hilbert(temp)));
title('\phi(r,t)');
hold on
%ylim([((150-1)*dz*1e3+16.5) ((180-1)*dz*1e3+16.5)]);
ylim([125 205]);
%ylim([125 205]);
set(gcf, 'Position', get(0, 'Screensize'));
plot(col1, row1, 'r.')
% ylabel(' Depth');
% xlabel(' Lateral');
%dat = data(end - ceil(13.5e-3/ dz):end- ceil(3.5e-3/ dz), Nx/2 - ceil(5e-3/dx) : Nx/2 +ceil(5e-3/dx) ,:);
figure;
imagesc(abs(hilbert(dat(:,:,i))));
set(gcf, 'Position', get(0, 'Screensize'));
hold on
plot(col1, row1, 'r.')
title('b(r,t)');
ylim([125 205]);
%ylim([150 180]);
%ylim([125 205]);
% ylabel(' Depth');
% xlabel(' Lateral');
%
%map = speed_map2(end - ceil(13.5e-3/ dz)-1:end- ceil(3.5e-3/ dz)-1, Nx/2 - ceil(5e-3/dx) : Nx/2 +ceil(5e-3/dx) ,:);
map = speed_map1(end - ceil(13.5e-3/ dz)-1:end- ceil(3.5e-3/ dz)-1, Nx/2 - ceil(5e-3/dx) : Nx/2 +ceil(5e-3/dx) ,:);
% gt1 = groundtruth1(end - ceil(13.5e-3/ dz)-1:end- ceil(3.5e-3/ dz)-1, Nx/2 - ceil(5e-3/dx) : Nx/2 +ceil(5e-3/dx),i);
% f = find(gt1);
% [I, J] = ind2sub(size(gt1),f);
figure;
imagesc(map);
set(gcf, 'Position', get(0, 'Screensize'));
% 
hold on
plot(col1, row1, 'r.')
title('Recovered');
%ylim([150 180]);
ylim([125 205]);
% ylabel(' Depth');
% xlabel(' Lateral');
colorbar
% figure;
% imagesc(map1);

% 
% hold on 
% plot(J,I,'k.')
% title('Groundtruth');
% %ylim([150 180]);
% set(gcf, 'Position', get(0, 'Screensize'));

%% part2
f0 = 5e6; % Transducer center frequency [Hz]
c = 1540; % Speed of sound [m/s]
lambda = c/f0;
dx = lambda/10;  
dz = lambda/10;
% figure;
% imagesc(zeros(326,327));
figure;
imagesc(map);
set(gcf, 'Position', get(0, 'Screensize'));
%imagesc(linspace(0,(a(2)-1)*dx*1e3,a(2)) -(a(2)-1)*dx*1e3/2,linspace(0,(a(1)-1)*dz*1e3 ,a(1))+ 19,zeros(164,327));

for i=30:120 
% temp = u(ceil(2.5e-3/ dz)+4-1:end -ceil(2.5e-3/ dz)+4 ,: ,i);
% temp = temp/max(max(temp));
% a =size(temp);
temp = u(: ,: ,i);
[~, row1, col1] = get_centroids_coef(temp);

%columnmodified = (col1-1) * dx*1e3 -5.0204  ;  
%rowmodified = (row1-1) * dz*1e3 +19 ;  

hold on
plot(col1, row1, 'k.')

end
ylim([125 205]);
%ylim([150 180]);
title('Superposed Image');
colorbar
% ylabel(' Depth');
% xlabel(' Lateral');
%% Cocentric 

% f0 = 5e6; % Transducer center frequency [Hz]
% c = 1540; % Speed of sound [m/s]
% lambda = c/f0; % Wave length [m]
% FrameRate =150; %Hz
% deltat = 1 / FrameRate;
% numberoftimesteps = 1*FrameRate +1; 
% dx = lambda/10; % 
% dz = lambda/10;
% sigma_b = lambda;
% Nz =975;
% Nx =1238;
% i = 75;
% temp = u(: ,: ,i);
% a =size(temp);
% figure;
% imagesc(abs(hilbert(temp)));
% xlim([350 600]);
% [~, row1, col1] = get_centroids_coef(temp);
% hold on
% plot(col1, row1, 'r.')
% 
% title('\phi(r,t)');
% figure;
% imagesc(abs(hilbert(dat(:,:,i))));
% 
% hold on
% plot(col1, row1, 'r.')
% title('b(r,t)');
% xlim([350 600]);
% map = angular(end - ceil(13.5e-3/ dz)-1:end - ceil(5e-3/ dz)-1, Nx/2 - ceil(10e-3/dx): Nx/2 +ceil(10e-3/dx));
% figure;
% imagesc(map);
% %set(gcf, 'Position', get(0, 'Screensize'));
% % 
% hold on
% plot(col1, row1, 'r.')
% xlim([350 600]);
% title('Recovered');
% %% part2
% 
% f0 = 5e6; % Transducer center frequency [Hz]
% c = 1540; % Speed of sound [m/s]
% lambda = c/f0;
% dx = lambda/10;  
% dz = lambda/10;
% % figure;
% % imagesc(zeros(326,327));
% figure;
% imagesc(map);
% set(gcf, 'Position', get(0, 'Screensize'));
% %imagesc(linspace(0,(a(2)-1)*dx*1e3,a(2)) -(a(2)-1)*dx*1e3/2,linspace(0,(a(1)-1)*dz*1e3 ,a(1))+ 19,zeros(164,327));
% 
% for i=40:90
% 
% temp = u(: ,: ,i);
% [~, row1, col1] = get_centroids_coef(temp);
% 
% hold on
% plot(col1, row1, 'r.')
% 
% end
% %ylim([125 205]);
% %ylim([150 180]);
% title('Superposed Image');
% colorbar
% xlim([350 600]);
% % ylabel(' Depth');
% % xlabel(' Lateral');
