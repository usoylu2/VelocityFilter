%% OppositeDirectionVessels 
% f0 = 15e6; % Transducer center frequency [Hz]
% c = 1540; % Speed of sound [m/s]
% lambda = c/f0; % Wave length [m]
% FrameRate =100; %Hz
% deltat = 1 / FrameRate;
% numberoftimesteps = 1*FrameRate +1; 
% dx = lambda/10; % 
% dz = lambda/10;
% sigma_b = lambda;
% Nz =975;
% Nx =1841;
% 
% startz = 0; % m
% endz = (Nz)*dz; %m
% startx = 0; % m
% endx = (Nx)*dx; %m
% kG =3/lambda;
% sigma_w =1; %seconds
% v0 = [1e-3,0]; %x and z
% 
% t = -1:deltat :1 ;
% filter = zeros(Nz,Nx,length(t));
% J = zeros(Nz,Nx,length(t));
% 
% for i=1:length(t)
%     filter(:,:,i) =exp(-1/2 *(t(i)/sigma_w^2)^2);
% end
% 
% z = linspace(dz/2,endz - dz/2,Nz);
% x = linspace(dx/2,endx - dx/2,Nx);
% [Linex,LineZ] = meshgrid(x,z);
% centerz = 9e-3;
% 
% for i=1:length(t)
%     J(:,:,i) =deltat * dx * dz*(kG /(sqrt(2 *pi)* sigma_w )) * besselj(1,2*pi*kG * sqrt((LineZ-centerz-v0(2)*t(i)).^2 + (Linex-endx/2 -v0(1)*t(i)).^2))./ (4 * pi^2 *  sqrt((LineZ-centerz -v0(2)*t(i)).^2 + (Linex-endx/2- v0(1)*t(i)).^2) );
% end
% 
% filter = filter.* J;
% %%
% %GT=groundtruth(end - ceil(2e-3/ dz)-1:end,:);
% %PSF = data(end - ceil(2e-3/ dz)-1:end,920-500:1:920+500);
% filt = filter(end - ceil(1.5e-3/ dz)-1:end- ceil(0.5e-3/ dz)-1, Nx/2 - ceil(5e-3/dx)-1 : Nx/2 +ceil(5e-3/dx)-1 ,:);
% dat = data(end - ceil(1.5e-3/ dz)-1:end- ceil(0.5e-3/ dz)-1, Nx/2 - ceil(5e-3/dx)-1 : Nx/2 +ceil(5e-3/dx)-1 ,51:151);
% filt =single(filt);
% dat =single(dat);
% u = convnfft(dat,filt,'same');

%% Intersecting Vessels
% 
% f0 = 15e6; % Transducer center frequency [Hz]
% c = 1540; % Speed of sound [m/s]
% lambda = c/f0; % Wave length [m]
% FrameRate =100; %Hz
% deltat = 1 / FrameRate;
% numberoftimesteps = 1*FrameRate +1; 
% dx = lambda/10; % 
% dz = lambda/10;
% sigma_b = lambda;
% Nz =975;
% Nx =1841;
% 
% startz = 0; % m
% endz = (Nz)*dz; %m
% startx = 0; % m
% endx = (Nx)*dx; %m
% kG =3/lambda;
% sigma_w =1; %seconds
% v= 9.9e-3;
% slope = -0.1;
% v0 = [v*cos(abs(atan(slope))),v*sin(abs(atan(slope)))]; %x and z
% 
% t = -1:deltat :1 ;
% filter = zeros(Nz,Nx,length(t));
% J = zeros(Nz,Nx,length(t));
% 
% for i=1:length(t)
%     filter(:,:,i) =exp(-1/2 *(t(i)/sigma_w^2)^2);
% end
% 
% z = linspace(dz/2,endz - dz/2,Nz);
% x = linspace(dx/2,endx - dx/2,Nx);
% [Linex,LineZ] = meshgrid(x,z);
% centerz = 8e-3;
% 
% for i=1:length(t)
%     J(:,:,i) =deltat * dx * dz*(kG /(sqrt(2 *pi)* sigma_w )) * besselj(1,2*pi*kG * sqrt((LineZ-centerz-v0(2)*t(i)).^2 + (Linex-endx/2 -v0(1)*t(i)).^2))./ (4 * pi^2 *  sqrt((LineZ-centerz -v0(2)*t(i)).^2 + (Linex-endx/2- v0(1)*t(i)).^2) );
% end
% 
% filter = filter.* J;
% %%
% %GT=groundtruth(end - ceil(2e-3/ dz)-1:end,:);
% %PSF = data(end - ceil(3e-3/ dz)-1:end-ceil(1e-3/ dz)-1,920-500:1:920+500);
% filt = filter(end - ceil(3.5e-3/ dz)-1:end- ceil(0.5e-3/ dz)-1, Nx/2 - ceil(5e-3/dx)-1 : Nx/2 +ceil(5e-3/dx)-1 ,:);
% dat = data(end - ceil(3.5e-3/ dz)-1:end- ceil(0.5e-3/ dz)-1, Nx/2 - ceil(5e-3/dx)-1 : Nx/2 +ceil(5e-3/dx)-1 ,51:151);
% C = zeros(293,977,300);
% C(:,:,100:200)=dat;
% C =single(C);
% filt =single(filt);
% dat =single(dat);
% u = convnfft(C,filt,'same');
% result = u(:,:,100:200);

%% Cocentric Vessels

f0 = 5e6; % Transducer center frequency [Hz]
c = 1540; % Speed of sound [m/s]
lambda = c/f0; % Wave length [m]
FrameRate =100; %Hz
deltat = 1 / FrameRate;
dx = lambda/10; % 
dz = lambda/10;
sigma_b = lambda;
width = 0.27e-3; % Width of element
element_height = 5/1000; % Height of element [m]
kerf = 0.03/1000; % Kerf [m]
chano = 128; 
startz = 0; % m
endz = 30e-3; %m
x = -(chano-1)/2*(width+kerf):dx:(chano-1)/2*(width+kerf);
z = startz:dz:endz; 
[Linex,LineZ] = meshgrid(x,z);
centerz = 21.5e-3;
Nz =length(z);
Nx =length(x);

kG =3/lambda;
sigma_w =0.4; %seconds
v= 10e-3;
angle = 7*pi/8;
v0 = [v*sin(angle),-v*cos(angle)]; %x and z
%v0 =[v,0];

t = -0.5:deltat :0.5 ;
filter = zeros(Nz,Nx,length(t));
J = zeros(Nz,Nx,length(t));

for i=1:length(t)
    filter(:,:,i) =exp(-1/2 *(t(i)/sigma_w^2)^2);
end

for i=1:length(t)
    J(:,:,i) = besselj(1,2*pi*kG * sqrt((LineZ-centerz-v0(2)*t(i)).^2 + (Linex -v0(1)*t(i)).^2))./ (  sqrt((LineZ-centerz -v0(2)*t(i)).^2 + (Linex - v0(1)*t(i)).^2) );
end

filter =  filter.* J;
filter = (kG /(sqrt(2 *pi)* sigma_w )) *deltat * dx * dz* filter /(4 * pi^2 );
%%
filt = filter(end - ceil(16e-3/ dz)-1:end- ceil(1e-3/ dz)-1 , Nx/2 - ceil(15e-3/dx)-1 : Nx/2 +ceil(15e-3/dx)-1 ,:);
% dat = data(end - ceil(16e-3/ dz)-1:end - ceil(1e-3/ dz)-1, Nx/2 - ceil(15e-3/dx)-1 : Nx/2 +ceil(15e-3/dx)-1 ,:);

filt =single(filt);
dat =single(dat);
filter = [];
x =[];
z =[];
Linex =[];
LineZ =[];
u = convnfft(dat,filt,'same');
save('filterat7piover8','u');