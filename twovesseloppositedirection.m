clear all;close all;clc;
% path(path, 'C:\Users\usoylu2\Desktop\Field_II_ver_3_24_windows')
% field_init
% set_field ('att',1.5*100); 
% set_field ('Freq_att',0.5*100/1e6);
% set_field ('att_f0',5e6); 
% set_field ('use_att',1);
%% Random microbubles & Field II simulation
f0 = 15e6; % Transducer center frequency [Hz]
fs = 100e6; % Sampling frequency [Hz]
c = 1540; % Speed of sound [m/s]
lambda = c/f0; % Wave length [m]
width = 0.27e-3; % Width of element
element_height = 5/1000; % Height of element [m]
kerf = 0.03/1000; % Kerf [m]
focus=[0 0 20]/1000; % Fixed focal point [m]
N_elements = 64; % Number of elements in the transducer
N_active = N_elements; % Active elements in the transducer
% Set the sampling frequency
set_sampling(fs);
% Generate aperture for emission
emit_aperture = xdc_linear_array(N_active, width, element_height, kerf, 10, 10, focus);
% Set the impulse response and excitation of the emit aperture
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response.*hann(max(size(impulse_response)))';
xdc_impulse(emit_aperture, impulse_response);
excitation = sin(2*pi*f0*(0:1/fs:1/f0));
excitation = zeros(length(excitation),1)';
excitation(9) = 1;
xdc_excitation(emit_aperture, excitation);
% Generate aperture for reception
% receive_aperture = xdc_convex_array (N_elements, width, element_height, kerf, Rconvex, 5, 5, focus);
receive_aperture = xdc_linear_array(N_active, width, element_height, kerf, 10, 10, focus);
% Set the impulse response for the receive aperture
xdc_impulse(receive_aperture, impulse_response);

FrameRate = 100; %Hz
duration = 2*FrameRate +1; %2 seconds
deltat = 1 / FrameRate;

chano = N_elements;
resol = 1*lambda;  %axial resolution
dx = resol/10; % lateral resolution
dz = resol/10;
startdepth = 0; % m
enddepth = 10e-3; %m

x = -(chano-1)/2*(width+kerf):dx:(chano-1)/2*(width+kerf); % lateral dimension of the FOV
z = startdepth:dz:enddepth; % axial dimension of the FOV
arrayx = (-chano/2+.5:chano/2-.5)*(width+kerf);
arrayz = zeros(1,length(arrayx));
senscutoff = 0.0;
Origin = [0,0];
Nt = 1;
[Linex,LineZ] = meshgrid(x,z);

Nz = length(z);
Nx = length(x);
data = zeros(Nz,Nx,duration);
centerlocation =9e-3;

%Parabolic Profile
D = 0.1e-3; %diameter of the vessels
z1 = centerlocation- D/2;
z2 = centerlocation + D/2;
v1 = 10e-3; %max velocites in millimeters/seconds
v2 = -10e-3;
parabol1 = (1 - (LineZ - z1).^2 /(D/2)^2);
parabol2 = (1 - (LineZ - z2).^2 /(D/2)^2);

%Number of Bubbles

density = 300; %bubbles per mm3
len = abs(x(1)-x(end))* 1e3; %mm

numberofbubbles = 2* pi * (D/2 *1e3)^2 * len * density;
numberofbubbles = ceil(numberofbubbles);
%Paralel Lines

slope = 0;
line1 = (((LineZ - z1)+ slope*(Linex ) >= (-D/2))...
    &((LineZ -  z1)+ slope*(Linex ) <= (D/2)));

line2 = (((LineZ - z2)+ slope*(Linex ) >= (-D/2))...
    &((LineZ -  z2)+ slope*(Linex ) <= (D/2)));

linea = line1 .* parabol1; 
lineb = line2 .* parabol2;
speed_map1 = ones(Nz,Nx);
speed_map1 = v1*speed_map1;
speed_map1 = speed_map1 .* linea;

speed_map2 = ones(Nz,Nx);
speed_map2 = v2*speed_map2;
speed_map2 = speed_map2 .* lineb;

speed_map = speed_map1 +speed_map2;
directionx = zeros(Nz,Nx);
directionz = zeros(Nz,Nx);
directionx(:,:) = sqrt(1/ (1+slope^2)); 
directionz(:,:) = directionx * slope;

lineorig = line1 +line2;
directionx = directionx .* lineorig;
directionz = directionz .* lineorig;
velocityx = directionx .* speed_map;
velocityz = -1*directionz .* speed_map;

positions = rand(numberofbubbles, 2);
positions(:,1) = positions(:,1)*abs(x(1)-x(end));
positions(:,1) = positions(:,1)-abs(x(1)-x(end))/2;
positions(:,2) = positions(:,2)*2*D;
positions(:,2) = positions(:,2) + centerlocation -D;

directionx = [];
directionz = [];
Linex =[];
LineZ = [];
linea = [];
lineb = [];
line1 = [];
line2 = [];
parabol1 = [];
parabol2 = [];
speed_map =[];
speed_map1 =[];
speed_map2 =[];

for p=numberofbubbles:-1:1
    if lineorig(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /(abs(x(1)-x(end))) * Nx))==0
        positions(p,:) = [];
    end
end

temp = size(positions);
numberofbubbles = temp(1);

for number_of_frames=1:duration
fprintf("Frame Number is %d",number_of_frames);

bublepositions = zeros(numberofbubbles,3);
bublepositions(:,1) = positions(:,1);
bublepositions(:,3) = positions(:,2);
phantom_positions = bublepositions;
phantom_amplitudes = ones(numberofbubbles,1)*1e25;

senscutoff = 0.0;
Origin = [0,0];
Nt = 1;
[X,Z] = meshgrid(x,z);
[Na,Nl] = size(X);

X = repmat(X,[1,1,chano]);
Z = repmat(Z,[1,1,chano]);
arrayX = repmat(reshape(arrayx,[1,1,chano]),[Na,Nl,1]);% Meshgrid the array vector
arrayZ = repmat(reshape(arrayz,[1,1,chano]),[Na,Nl,1]);% Meshgrid the array vector
twpeak = 0.20e-6; %time offset for true ultrasound peak arrival time.

scat =[];
[scat,start_time] = calc_scat_all(emit_aperture,receive_aperture,phantom_positions,phantom_amplitudes,1);
sum(sum(sum(isnan(scat))))~= 0
scat = [zeros(round(start_time*fs),N_active*N_active);scat];

%%
time_samples = size(scat);
scat = reshape(scat,[time_samples(1),N_elements,N_elements]);
scat = sum(scat,2);
scat = squeeze(scat);
% scat = scat';
%% DAS beamforming
if sum(sum(sum(isnan(scat))))~= 0
    disp('NAN. Why?');
else
    disp('No problem.'); 
    Steer = 0 /180 *pi;
    NoAngles = length(Steer);
    N_transmit = N_elements;
    N_receive = N_elements;
    RData = scat;
    chano = N_receive;
    Fs = fs; %sampling frequency
    ft = f0; %transmit frequency
    cycle = 2; % number of transmit cycles
    plength = 1/ft*cycle; % transmit pulse length (sec)
    plengths = ceil(plength*Fs); % transmit pulse length in unit of samples
    % twpeak = 3.5781/fc;
    % prf = 500;

    PageEnd = size(RData,1);
    % Reshape the raw channel data to 3D (x,z,t)
    ReShpRData = zeros(PageEnd,chano,NoAngles,Nt);
    Tx_focus = 0e-3; %m
    if Tx_focus > 0 % not plane waves
        FocalPt(1) = Origin(1) + Tx_focus * sin(Steer);
        FocalPt(2) = 0.0;
        FocalPt(3) = Tx_focus * cos(Steer);
        % Compute distance to focal point from each active element.
        X_array = arrayx - FocalPt(1);
        Tx_Delay = sqrt(X_array.*X_array + FocalPt(3)*FocalPt(3))/c;
        Tx_Delay = max(Tx_Delay(:)) - Tx_Delay;
    elseif Tx_focus == 0   %plane wave imaging
        Tx_Delay = (arrayx'- Origin(1)) * sin(Steer)/c;
        Tx_Delay = Tx_Delay - min(Tx_Delay(:));
    else   %diverging waves
        FocalPt(1) = Origin(1) + Tx_focus * sin(Steer);
        FocalPt(2) = 0;
        FocalPt(3) = Tx_focus * cos(Steer);
        % Compute distance to focal point from each active element.
        X_array = arrayx' - FocalPt(1);
        Tx_Delay = sqrt(X_array.*X_array + FocalPt(3)*FocalPt(3))/c;
        Tx_Delay = Tx_Delay - min(Tx_Delay(:));
    end
    for i = 1:Nt
        for j = 1:NoAngles
            dummy = double(squeeze(RData((i-1)*NoAngles*PageEnd+(j-1)*PageEnd+1:(i-1)*NoAngles*PageEnd+j*PageEnd,:)));
            ReShpRData(:,:,j,i) = reshape(dummy,PageEnd,N_receive,1,1);
        end
    end
    %ReShpRData = scat for single plane wave imaging

    % Pixel-oriented beamforming (matrix version)

    clear delaysub delayidx;
    for n = 1:NoAngles
        Tdelay = interp1(arrayx,Tx_Delay,x,'linear');
        Tdelay = repmat(reshape(Tdelay,[1,Nl,1]),[Na,1,chano]); % Transmit delay
%         delay = round(((Z+sqrt((Z-arrayZ).^2+(X-arrayX).^2))/c+Tdelay+lensCorrection/c*2+twpeak)*Fs);% calculate the delay in each channel for each pixel in units of samples
        delay = round(((Z+sqrt((Z-arrayZ).^2+(X-arrayX).^2))/c+Tdelay+twpeak)*Fs);
        delay(delay<=0) = 1; % clear 0 delays
        delay(delay>PageEnd) = PageEnd; % clear delays beyond the sampling range
        %delaysub(:,:,:,n) = delay;
        %delayidx(:,:,:,n) = sub2ind([size(ReShpRData,1),chano],delay,repmat(reshape(1:chano,[1,1,chano]),[Na,Nl,1]));%find the indices of the data for each pixel in each channel
        delay = sub2ind([size(ReShpRData,1),chano],delay,repmat(reshape(1:chano,[1,1,chano]),[Na,Nl,1]));%find the indices of the data for each pixel in each channel
    end
        
    X =[];
    Z =[];
    arrayZ =[];
    arrayX =[];
    Tdelay =[];
    % Beamforming and IQ downmixing
    beamformedIQ = zeros(Na,Nl,NoAngles,Nt);

    for i = 1:Nt
        for n = 1:NoAngles
            dummy = squeeze(ReShpRData(:,:,n,i));
            %size(dummy) is 6547*128
            beamformedIQ = sum(dummy(delay),3);
            %TEMP = sum(dummy(delayidx(:,:,:,n)).*ElementSens,3);
            %size(TEMP) is 1624*495, size(delayidx is 1624*495*128)
            %beamformedIQ(:,:,n,i) = TEMP;
            %env(:,:,n,i) = abs(hilbert(TEMP));
        end
    end
    %env = beamformedIQ;
    %data(:,:) = beamformedIQ;
    %figure;imagesc(x*1e3,z*1e3,20*log10(env./max(env(:))),[-60 0]);colormap(gray);
    %figure;imagesc(x*1e3,z*1e3,beamformedIQ);colormap(gray);

    beamformedIQ = abs( hilbert(beamformedIQ));
    %figure;imagesc(x*1e3,z*1e3,beamformedIQ);colormap(gray);
    data(:,:,number_of_frames) =beamformedIQ;
    %save(strcat('aupone',int2str(u),'.mat'),'beamformedIQ');
    %save(strcat('auponegt',int2str(u),'.mat'),'groundtruth');
    
    %Update positions

    for p=1:numberofbubbles
        delvx = velocityx(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /(abs(x(1)-x(end))) * Nx));
        delvz = velocityz(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /(abs(x(1)-x(end))) * Nx));
        positions(p,:) = [delvx,delvz]*deltat + positions(p,:);
    end
    

%Delete the bubbles which are out of the scene 
for p=numberofbubbles:-1:1
    if positions(p,1) < -abs(x(1)-x(end))/2 || positions(p,1) > abs(x(1)-x(end))/2 || positions(p,2) < 2e-3 || positions(p,2) > enddepth
        positions(p,:) = [];
    end
end
temp = size(positions);
numberofbubbles = temp(1);
%Delete the bubbles which are not moving
for p=numberofbubbles:-1:1
    if lineorig(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /(abs(x(1)-x(end))) * Nx))==0
        positions(p,:) = [];
    end
end
temp = size(positions);
numberofbubbles = temp(1);

end
end
%save('twovesseloppositedirection.mat','data','-v7.3');
%% 
% f0 = 5e6; % Transducer center frequency [Hz]
% c = 1540; % Speed of sound [m/s]
% lambda = c/f0; % Wave length [m]
% FrameRate =100; %Hz
% deltat = 1 / FrameRate;
% numberoftimesteps = 2*FrameRate +1; 
% dx = lambda/4; % 
% dz = lambda/10;
% sigma_b = lambda;
% Nz =1000;
% Nx = 400;
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
% filter = zeros(Nz,Nx,len(t));
% J = zeros(Nz,Nx,len(t));
% 
% for i=1:len(t)
%     filter(:,:,i) =exp(-1/2 *(t(i)/sigma_w^2)^2);
% end
% 
% z = linspace(dz/2,endz - dz/2,Nz);
% x = linspace(dx/2,endx - dx/2,Nx);
% [Linex,LineZ] = meshgrid(x,z);
% 
% 
% for i=1:len(t)
%     J(:,:,i) =deltat * dx * dz*(kG /(sqrt(2 *pi)* sigma_w )) * besselj(1,2*pi*kG * sqrt((LineZ-endz/2-v0(2)*t(i)).^2 + (Linex-endx/2 -v0(1)*t(i)).^2))./ (4 * pi^2 *  sqrt((LineZ-endz/2 -v0(2)*t(i)).^2 + (Linex-endx/2- v0(1)*t(i)).^2) );
% end
% 
% filter = filter.* J;
% 
% %%
% 
% filt = filter(470:530,170:230,:);
% dat = data(712:912,125:end-125 ,:);
% u = convnfft(dat,filt,'same');