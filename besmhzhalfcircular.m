clear all;close all;clc;
% path(path, 'C:\Users\usoylu2\Desktop\Field_II_ver_3_24_windows')
% field_init
% set_field ('att',1.5*100); 
% set_field ('Freq_att',0.5*100/1e6);
% set_field ('att_f0',5e6); 
% set_field ('use_att',1);
%% Random microbubles & Field II simulation
f0 = 5e6; % Transducer center frequency [Hz]
fs = 100e6; % Sampling frequency [Hz]
c = 1540; % Speed of sound [m/s]
lambda = c/f0; % Wave length [m]
width = 0.27e-3; % Width of element
element_height = 5/1000; % Height of element [m]
kerf = 0.03/1000; % Kerf [m]
focus=[0 0 20]/1000; % Fixed focal point [m]
N_elements = 128; % Number of elements in the transducer
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
duration = 1*FrameRate+1; %2 seconds
deltat = 1/FrameRate;

chano = N_elements;
resol = 1*lambda;  %axial resolution
dx = resol/10; % lateral resolution
dz = resol/10;
startdepth = 0; % m
enddepth = 30e-3; %m

x = -(chano-1)/2*(width+kerf):dx:(chano-1)/2*(width+kerf); % lateral dimension of the FOV
z = startdepth:dz:enddepth; % axial dimension of the FOV
arrayx = (-chano/2+.5:chano/2-.5)*(width+kerf);
arrayz = zeros(1,length(arrayx));
senscutoff = 0.0;
Origin = [0,0];
Nt = 1;
[Linex,LineZ]=meshgrid(x,z);

numberofvessels =3;
Nz =length(z);
Nx = length(x);
datatemp = zeros(Nz,Nx,duration,numberofvessels);

D = 0.3e-3;
% r1 = 14e-3;
% r2 = 10e-3;
% r3 = 6e-3;
r1 = 12e-3;
r2 = 11.7e-3;
r3 = 6e-3;
centerz1 = 29e-3;
centerz2 =centerz1 ;
centerz3 = 29e-3;
%Circular Bloood Structure
circle1 = (((LineZ-centerz1).^2+(Linex).^2 >= (r1 - D/2)^2)...
    &((LineZ-centerz1).^2+(Linex).^2 <= (r1 + D/2)^2));
parabol1 = (1 -abs(sqrt((LineZ-centerz1).^2+(Linex).^2) - r1 ).^2 /(D/2)^2);

circle2 = (((LineZ-centerz2).^2+(Linex).^2 >= (r2 - D/2)^2)...
    &((LineZ-centerz2).^2+(Linex).^2 <= (r2 + D/2)^2));
parabol2 = (1 -abs(sqrt((LineZ-centerz2).^2+(Linex).^2) - r2 ).^2 /(D/2)^2);

circle3 =(((LineZ-centerz3).^2+(Linex).^2 >= (r3 - D/2)^2)...
   &((LineZ-centerz3).^2+(Linex).^2 <= (r3 + D/2)^2));
parabol3 = (1 -abs(sqrt((LineZ-centerz3).^2+(Linex).^2) - r3 ).^2 /(D/2)^2);

speed = 10e-3;
angular_velocity1 = speed./sqrt((LineZ-centerz1).^2+(Linex).^2); 
angular_velocity1 =angular_velocity1 .*circle1.*parabol1;
theta1 = angular_velocity1 .* deltat;

speed = 10e-3;
angular_velocity2 = speed./sqrt((LineZ-centerz2).^2+(Linex).^2); 
angular_velocity2 =angular_velocity2 .*circle2.*parabol2*(-1);
theta2 = angular_velocity2 .* deltat ;

speed = 10e-3;
angular_velocity3 = speed./sqrt((LineZ-centerz3).^2+(Linex).^2); 
angular_velocity3 =angular_velocity3 .*circle3.*parabol3;
theta3 = angular_velocity3 .* deltat;

len1 = pi * r1 *1e3;
len2 = pi * r2 *1e3;
len3 = pi * r3 *1e3;

numberofbubbles1 =  pi * (D/2 *1e3)^2 * len1 * 100;
numberofbubbles1 = ceil(numberofbubbles1) ;
xu = rand(numberofbubbles1, 1);
xu = xu*(2*r1 + D) - (r1 +D/2);
d = rand(numberofbubbles1, 1);
d = d * (D) + (r1 - D/2);
zu = -1*sqrt(abs(d.^2 - xu.^2)) +centerz1;
positions1(:,1)=xu;
positions1(:,2)=zu;

for p=numberofbubbles1:-1:1
    if positions1(p,1) < -(abs(x(1)-x(end)))/2|| positions1(p,1) > (abs(x(1)-x(end)))/2|| positions1(p,2) < 2e-3 || positions1(p,2) > (enddepth-0.1e-3)
        positions1(p,:) = [];
    end
end
temp = size(positions1);
numberofbubbles1 = temp(1);    

for p=numberofbubbles1:-1:1
    disp(p);
    if circle1(ceil(positions1(p,2)/enddepth * Nz),ceil((positions1(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx))==0
        positions1(p,:) = [];
    end
end
temp = size(positions1);
numberofbubbles1 = temp(1);

numberofbubbles2 =  pi * (D/2 *1e3)^2 * len2 * 100;
numberofbubbles2 = ceil(numberofbubbles2) ;
xr = rand(numberofbubbles2, 1);
xr = xr*(2*r2 + D) - (r2 +D/2);
dr = rand(numberofbubbles2, 1);
dr = dr * (D) + (r2 - D/2);
zr =-1*sqrt(abs(dr.^2 - xr.^2)) +centerz2;
positions2(:,1)=xr;
positions2(:,2)=zr;

for p=numberofbubbles2:-1:1
    if positions2(p,1) < -(abs(x(1)-x(end)))/2|| positions2(p,1) > (abs(x(1)-x(end)))/2|| positions2(p,2) < 2e-3 || positions2(p,2) > (enddepth-0.1e-3)
        positions2(p,:) = [];
    end
end
temp = size(positions2);
numberofbubbles2 = temp(1);  

for p=numberofbubbles2:-1:1
    disp(p);
    if circle2(ceil(positions2(p,2)/enddepth * Nz),ceil((positions2(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx))==0
        positions2(p,:) = [];
    end
end
temp = size(positions2);
numberofbubbles2 = temp(1);

numberofbubbles3 =  pi * (D/2 *1e3)^2 * len3 * 100;
numberofbubbles3 = ceil(numberofbubbles3) ;
xg = rand(numberofbubbles3, 1);
xg = xg*(2*r3 + D) - (r3 +D/2);
gr = rand(numberofbubbles3, 1);
gr = gr * (D) + (r3 - D/2);
zg =-1*sqrt(abs(gr.^2 - xg.^2)) +centerz3;
positions3(:,1)=xg;
positions3(:,2)=zg;

for p=numberofbubbles3:-1:1
    if positions3(p,1) < -(abs(x(1)-x(end)))/2|| positions3(p,1) > (abs(x(1)-x(end)))/2|| positions3(p,2) < 2e-3 || positions3(p,2) > (enddepth-0.1e-3)
        positions3(p,:) = [];
    end
end
temp = size(positions3);
numberofbubbles3 = temp(1);  

for p=numberofbubbles3:-1:1
    disp(p);
    if circle3(ceil(positions3(p,2)/enddepth * Nz),ceil((positions3(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx))==0
        positions3(p,:) = [];
    end
end
temp = size(positions3);
numberofbubbles3 = temp(1);

angular_velocity1 = [];
angular_velocity2 = [];
angular_velocity3 = [];
Linex =[];
LineZ = [];
parabol3 = [];
parabol1 = [];
parabol2 = [];

for number_of_frames=1:duration
tic
fprintf("Frame Number is %d",number_of_frames); 

    for vessels=1:numberofvessels
     
        if vessels ==1
           positions = positions1; 
           theta =theta1;
           circle = circle1;
           numberofbubbles = numberofbubbles1;
           centerz = centerz1;
        elseif vessels ==2
           positions =positions2;
           theta =theta2;
           circle = circle2;
           numberofbubbles = numberofbubbles2;
           centerz = centerz2;
        elseif vessels ==3
           positions =positions3;
           theta =theta3;
           circle = circle3;
           numberofbubbles = numberofbubbles3;
           centerz = centerz3;
        end

bublepositions = zeros(numberofbubbles,3);
bublepositions(:,1) = positions(:,1);
bublepositions(:,3) = positions(:,2);
phantom_positions = bublepositions;
phantom_amplitudes = ones(numberofbubbles,1)*1e25;


arrayx = (-chano/2+.5:chano/2-.5)*(width+kerf);
arrayz = zeros(1,length(arrayx));
senscutoff = 0.0;
Origin = [0,0];
Nt = 1;
[X,Z] = meshgrid(x,z);
[Na,Nl] = size(X);

X = repmat(X,[1,1,chano]);
Z = repmat(Z,[1,1,chano]);
arrayX = repmat(reshape(arrayx,[1,1,chano]),[Na,Nl,1]);% Meshgrid the array vector
arrayZ = repmat(reshape(arrayz,[1,1,chano]),[Na,Nl,1]);% Meshgrid the array vector
twpeak = 0.50e-6; %time offset for true ultrasound peak arrival time.

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

    %beamformedIQ = abs( hilbert(beamformedIQ));
    %figure;imagesc(x*1e3,z*1e3,beamformedIQ);colormap(gray);
    datatemp(:,:,number_of_frames,vessels) =beamformedIQ;
    %save(strcat('aupone',int2str(u),'.mat'),'beamformedIQ');
    
    %Update positions

    for p=1:numberofbubbles
        r = sqrt((positions(p,2)-centerz).^2+(positions(p,1)).^2);
        if positions(p,2)-centerz <= 0 && positions(p,1)<=0
            beta = atan(abs(positions(p,2)-centerz) / abs(positions(p,1)));
            alpha = pi/2 -(beta+ theta(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx)));
            positions(p,1) = -r * sin(alpha);
            positions(p,2) = -r * cos(alpha) +centerz;
        elseif positions(p,2)-centerz< 0 && positions(p,1)>0
            beta = atan( abs(positions(p,1)) /abs(positions(p,2)-centerz));
            beta = beta + pi/2;
            alpha = pi/2 -(beta+ theta(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx)));
            positions(p,1) = -r * sin(alpha);
            positions(p,2) = -r * cos(alpha) +centerz;
        elseif positions(p,2)-centerz >=0 && positions(p,1)>=0
            beta = atan(abs(positions(p,2)-centerz) / abs(positions(p,1)));
            beta = beta + pi;
            alpha = pi/2 -(beta+ theta(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx)));
            positions(p,1) = -r * sin(alpha);
            positions(p,2) = -r * cos(alpha) +centerz;
        elseif positions(p,2)-centerz >0 && positions(p,1)<0
            beta = atan(abs(positions(p,1)) /abs(positions(p,2)-centerz));
            beta = beta + 3*pi/2;
            alpha = pi/2 -(beta+ theta(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx)));
            positions(p,1) = -r * sin(alpha);
            positions(p,2) = -r * cos(alpha) +centerz;
        end
    end
    
    
for p=numberofbubbles:-1:1
    if positions(p,1) < -(abs(x(1)-x(end)))/2|| positions(p,1) > (abs(x(1)-x(end)))/2|| positions(p,2) < 2e-3 || positions(p,2) > (centerz-0.1e-3)
        positions(p,:) = [];
    end
end
temp = size(positions);
numberofbubbles = temp(1);    
 
for p=numberofbubbles:-1:1
    disp(p);
    if circle(ceil(positions(p,2)/enddepth * Nz),ceil((positions(p,1)+abs(x(1)-x(end))/2) /abs(x(1)-x(end)) * Nx))==0
        disp('burda');
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
        elseif vessels ==3
           positions3 =positions;
           numberofbubbles3 = numberofbubbles;
        end       
end
end
toc
end
data = sum(datatemp,4);
save('halfcircle5MHzoppositerf','data','-v7.3');