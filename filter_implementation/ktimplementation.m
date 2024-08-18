%load('intersectingvessels5MHZrfwithgt150Hz.mat');
%Bubble Singal Information
f0 = 5e6; % Transducer center frequency [Hz]
c = 1540; % Speed of sound [m/s]
lambda = c/f0;
Dx = lambda /10;
Dz = lambda /10;
frame_rate = 150;% Hz
deltat = 1 / frame_rate;
Nx =327; %numbeofpixels
Nz =650; %numbeofpixels

t = -1:deltat :1 ;
kz = 2*pi*linspace(-1/(2*Dz),1/(2*Dz),Nz);
kx = 2*pi*linspace(-1/(2*Dx),1/(2*Dx),Nx);
[LineKx,LineKz,T] = meshgrid(kx,kz,t);
v0 = [0e-3,5e-3]; %x component and z component of velocity in meters per seconds
sigma_w = 0.5; 

filter = exp(-1i* (LineKx *v0(1) +LineKz*v0(2) ).*T ).* exp(-1/2 *(T/sigma_w^2).^2); %kz,kx,t
kx = [];
kz = [];
LineKx =[];
LineKz =[];
t = [];
T = [];
filter = single(filter);
%bubble = rand(Nz,Nx,450,'single');
%load('examplebubblefft2.mat');
load('shiftedbubblefft','bubble');
%bubble(:,:,451) = [];
tic
for i =1:Nx
    for j =1:Nz
        bubble(j,i,:)= fftmemoryonedim(bubble(j,i,:),filter(j,i,:),'same');
        %bubble(j,i,:)= overlappadd1D(bubble(j,i,:), filter(j,i,:),50);
    end
end
filter = [];
toc
bubble = ifftshift(bubble,1);
bubble = ifftshift(bubble,2);
bubble = ifft2(bubble);
bubble = real(bubble);
figure;
imagesc(abs(hilbert(bubble(:,:,225))));
colorbar;
%figure;
%imagesc(dat(:,:,225));
%colorbar;
