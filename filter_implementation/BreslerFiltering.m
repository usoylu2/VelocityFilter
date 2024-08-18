f0 = 5e6; % Transducer center frequency [Hz]
c = 1540; % Speed of sound [m/s]
lambda = c/f0; % Wave length [m]
FrameRate =150; %Hz
deltat = 1 / FrameRate;
numberoftimesteps = 1*FrameRate +1; 
dx = lambda/10; % 
dz = lambda/10;
sigma_b = lambda;
Nz =975;
Nx =1238;

startz = 0; % m
endz = (Nz)*dz; %m
startx = 0; % m
endx = (Nx)*dx; %m
kG =4*3/lambda; %or 4/lambda
sigma_w =0.5; %seconds
v =3e-3;
v0 = [v,0]; %x and z

t = -1:deltat :1 ;
filter = zeros(Nz,Nx,length(t));
J = zeros(Nz,Nx,length(t));

for i=1:length(t)
    filter(:,:,i) =exp(-1/2 *(t(i)/sigma_w^2)^2);
end

z = linspace(dz/2,endz - dz/2,Nz);
x = linspace(dx/2,endx - dx/2,Nx);
[Linex,LineZ] = meshgrid(x,z);
centerz = 21.5e-3;

for i=1:length(t)
    J(:,:,i) = deltat * dx * dz*(kG /(sqrt(2 *pi)* sigma_w )) * besselj(1,2*pi*kG * sqrt((LineZ-centerz-v0(2)*t(i)).^2 + (Linex-endx/2 -v0(1)*t(i)).^2))./ (4 * pi^2 *  sqrt((LineZ-centerz -v0(2)*t(i)).^2 + (Linex-endx/2- v0(1)*t(i)).^2));
end

filter = filter.* J;
filt = filter(end - ceil(13.5e-3/ dz):end- ceil(3.5e-3/ dz), Nx/2 - ceil(5e-3/dx) : Nx/2 +ceil(5e-3/dx) ,:);
%dat =data(end - ceil(13.5e-3/ dz):end- ceil(3.5e-3/ dz), Nx/2 - ceil(5e-3/dx): Nx/2 +ceil(5e-3/dx),:);
%dat =data(end - ceil(25e-3/ dz):end - ceil(5e-3/ dz), Nx/2 - ceil(7e-3/dx): Nx/2 +ceil(3e-3/dx),:);

filt =single(filt);
dat =single(dat);
filter =[];
x =[];
z =[];
Linex =[];
LineZ =[];
J =[];
t =[];
u = fftmemory(dat,filt,'same');
