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

FrameRate =100; %Hz
duration = 2*FrameRate+1; %2 seconds
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

Nz =length(z);
Nx = length(x);
data = zeros(Nz,Nx,duration);

%Paralel Lines
numberofvessels = 2;
datatemp = zeros(Nz,Nx,duration,numberofvessels);

%Parabolic Profile
centerline = 8e-3;
v1 = 10e-3; %max velocites in millimeters/seconds
D = 0.1e-3; %diameter of the vessels
v2 = 10e-3;
slope1 = 0;
slope2 = -0.2;
theta = atan(abs(slope2));
Dx = D /cos(theta);
line1 = (((LineZ - centerline)+ slope1*(Linex ) >= -(D/2))...
    &((LineZ - centerline)+ slope1*(Linex ) <= (D/2)));
line2 = (((LineZ - centerline)+ slope2*(Linex) >= -(Dx/2))...
    &((LineZ - centerline)+ slope2*(Linex ) <= (Dx/2)));

parabol1 = (1 - (LineZ - centerline).^2 /(D/2)^2);
parabol2 = (1 - (abs(LineZ- centerline+ slope2* Linex)/(sqrt(1 + slope2^2))).^2 /(D/2)^2);

linea = line1.*parabol1;
lineb = line2.*parabol2 ;
%Number of Bubbles

speed_map1 = ones(Nz,Nx);
directionx1 = zeros(Nz,Nx);
directionz1 = zeros(Nz,Nx);
directionx1(:,:) = sqrt(1/ (1+slope1^2)); 
directionz1(:,:) = directionx1 * slope1;
directionx1 =directionx1 .* line1 ;
directionz1 =directionz1 .* line1 ;

speed_map2 = ones(Nz,Nx);
directionx2 = zeros(Nz,Nx);
directionz2 = zeros(Nz,Nx);
directionx2(:,:) = sqrt(1/ (1+slope2^2)); 
directionz2(:,:) = directionx2 * slope2;
directionx2 =directionx2 .* line2 ;
directionz2 =directionz2 .* line2 ;

speed_map1 =v1*speed_map1;
speed_map1 = speed_map1 .* linea ;
velocityx1 = directionx1 .*speed_map1;
velocityz1 = -1*directionz1 .*speed_map1;

speed_map2 = v2*speed_map2;
speed_map2 = speed_map2 .* lineb ;
velocityx2 = directionx2 .*speed_map2;
velocityz2 = -1*directionz2 .*speed_map2;

velocityx = zeros(Nz,Nx);
velocityz = zeros(Nz,Nx);

len = (abs(x(1)-x(end)))*1e3; %mm
density = 300; %bubbles per mm3

numberofbubbles =  pi * (D/2 *1e3)^2 * len * density;
numberofbubbles = ceil(numberofbubbles);

positions1x = rand(numberofbubbles, 2);
positions1x(:,1) = positions1x(:,1)*(abs(x(1)-x(end)));
positions1x(:,1) = positions1x(:,1)-(abs(x(1)-x(end)))/2;
positions1x(:,2) = positions1x(:,2)*D+centerline -D/2;
positions1(:,2) = positions1x(:,2);
positions1(:,1) = positions1x(:,1);

for p=numberofbubbles:-1:1
    if line1(ceil(positions1(:,2)/enddepth * Nz),ceil((positions1(p,1)+(abs(x(1)-x(end)))/2) /(abs(x(1)-x(end))) * Nx))==0
        positions1(p,:) = [];
    end
end
temp = size(positions1);
numberofbubbles1 = temp(1);

len =sqrt( ((abs(x(1)-x(end))*1e3)* abs(slope2))^2+ (abs(x(1)-x(end))*1e3)^2  );
numberofbubbles =  pi * (D/2 *1e3)^2 * len * density ; %initilaze with more bubbles
numberofbubbles = ceil(numberofbubbles);
positionsx = rand(numberofbubbles, 1)*abs(x(1)-x(end))-abs(x(1)-x(end))/2;
d = rand(numberofbubbles,1)*Dx - Dx*0.5;
zu = -1 *slope2 * positionsx + d +centerline;

positions2(:,1) = positionsx;
positions2(:,2) = zu;

for p=numberofbubbles:-1:1
    if line2(ceil(positions2(p,2)/enddepth * Nz),ceil((positions2(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx))==0
        positions2(p,:) = [];
    end
end
temp = size(positions2);
numberofbubbles2 = temp(1);

directionx1 = [];
directionz1 = [];
directionx2 = [];
directionz2 = [];
Linex =[];
LineZ = [];
linea = [];
lineb = [];
parabol1 = [];
parabol2 = [];
speed_map =[];
speed_map1 =[];
speed_map2 =[];
positions1x =[];

for number_of_frames=1:duration
    fprintf("Frame Number is %d",number_of_frames);    
    for vessels=1:numberofvessels
     
        if vessels ==1
           positions = positions1; 
           velocityx =velocityx1;
           velocityz =velocityz1;
           line = line1;
           numberofbubbles = numberofbubbles1;
        elseif vessels ==2
           positions =positions2;
           velocityx =velocityx2;
           velocityz =velocityz2;
           line = line2;
           numberofbubbles = numberofbubbles2;
        end
        
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
%               delay = round(((Z+sqrt((Z-arrayZ).^2+(X-arrayX).^2))/c+Tdelay+lensCorrection/c*2+twpeak)*Fs);% calculate the delay in each channel for each pixel in units of samples
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
            datatemp(:,:,number_of_frames,vessels) = beamformedIQ;
            %save(strcat('aupone',int2str(u),'.mat'),'beamformedIQ');
            %save(strcat('auponegt',int2str(u),'.mat'),'groundtruth');
    
            %Update positions

            for p=1:numberofbubbles
                delvx = velocityx(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+(abs(x(1)-x(end)))/2) /(abs(x(1)-x(end))) * Nx));
                delvz = velocityz(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+(abs(x(1)-x(end)))/2) /(abs(x(1)-x(end))) * Nx));
                positions(p,:) = [delvx,delvz]*deltat + positions(p,:);
            end
    
        for p=numberofbubbles:-1:1
            if positions(p,1) < -(abs(x(1)-x(end)))/2|| positions(p,1) > (abs(x(1)-x(end)))/2|| positions(p,2) < 2e-3 || positions(p,2) > enddepth
                positions(p,:) = [];
            end
        end
        temp = size(positions);
        numberofbubbles = temp(1);
        %Delete the bubbles which are not moving
        for p=numberofbubbles:-1:1
            if line(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+(abs(x(1)-x(end)))/2) /(abs(x(1)-x(end))) * Nx))==0
                positions(p,:) = [];
            end
        end
        temp = size(positions);
        numberofbubbles = temp(1);

        if vessels ==1
           positions1 =positions; 
           numberofbubbles1 = numberofbubbles;
        elseif vessels ==2
           positions2 =positions;
           numberofbubbles2 = numberofbubbles;
        end        
        end
    end
end
data = sum(datatemp,4);
save('intersectingvessels.mat','data','-v7.3');
%save('datalinearintersectgt.mat','datagt');

