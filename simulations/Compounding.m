clear all;close all;clc;
%path(path,'/home/ufuk/Desktop/UltrasoundDL/Field_II_ver_3_22_linux/');
%field_init
% set_field ('att',1.5*100); 
% set_field ('Freq_att',0.5*100/1e6);
% set_field ('att_f0',5e6); 
% set_field ('use_att',1);
%% Random microbubles & Field II simulation
for u =1:1
tic
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
emit_aperture = xdc_linear_array (N_active, width, element_height, kerf, 10, 10, focus);
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

% Load the computer phantom
%phantom_positions = [0/1000,0/1000,10/1000;0/1000,0/1000,20/1000;0/1000,0/1000,30/1000;0/1000,0/1000,40/1000;...
%     -10/1000,0/1000,10/1000;-10/1000,0/1000,20/1000;-10/1000,0/1000,30/1000;-10/1000,0/1000,40/1000;...
%     10/1000,0/1000,10/1000;10/1000,0/1000,20/1000;10/1000,0/1000,30/1000;10/1000,0/1000,40/1000;];
%phantom_amplitudes = [1;1;1;1;1;1;1;1;1;1;1;1]*1e25;
%phantom_positions = [0/1000,0/1000,42/1000];
%phantom_amplitudes = [1]*1e25;
area = 19; %cm2 38mm * 50mm/100
density = 1;
number_of_bubles = area*density;
positions = rand(number_of_bubles, 2);
positions(:,1)= positions(:,1)*38e-3;
positions(:,1) = positions(:,1)-19e-3;
positions(:,2)= positions(:,2)*50e-3;
Nz =1624;
Nx =495;
groundtruth = zeros(Nz,Nx,'int8');
for j=1:number_of_bubles
    groundtruth(ceil(positions(j,2)/50e-3 * Nz ),ceil((positions(j,1)+19e-3) /38e-3 * Nx))=  groundtruth(ceil(positions(j,2)/50e-3 * Nz ),ceil((positions(j,1)+19e-3) /38e-3 * Nx))+1;
end
%data_gt(:,:,g+(u-1)*perdensity) = groundtruth; 
bublepositions = zeros(number_of_bubles,3);
bublepositions(:,1) = positions(:,1);
bublepositions(:,3) = positions(:,2);
phantom_positions = bublepositions;
phantom_amplitudes = ones(number_of_bubles,1)*1e25;

chano = N_elements;
resol = 1*lambda;  %axial resolution
dx = resol/4; % lateral resolution
dz = resol/10;
startdepth = 0; % m
enddepth = 50e-3; %m
x = -(chano-1)/2*(width+kerf):dx:(chano-1)/2*(width+kerf); % lateral dimension of the FOV
z = startdepth:dz:enddepth; % axial dimension of the FOV

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
%Theta = abs(atan((X-arrayX)./(Z-arrayZ))); % calculate the angle between each pixel to each element
%XT = width*lambda*pi*sin(Theta);
%ElementSens = abs(cos(Theta).*(sin(XT)./XT));% calculate element sensitivity (see VSX sequence programming)
% ElementSens(ElementSens<senscutoff) = 0; %set sensitivity cutoff
% ElementSens(ElementSens>=senscutoff) = 1;
twpeak = 0.50e-6; %time offset for true ultrasound peak arrival time.
% [~,cchanel] = min(abs(repmat(x,[chano,1]) - repmat(arrayx',[1,Nl])));%looking for closest element for each imaging line
data = zeros(Na,Nl);

% % Set the active elements using the apodization;
% apo = ones(1,N_active);
% apo=hann(N_active)';
%Tx_apo = hann(N_active)'; %Transmit apodization
%Rcv_apo = ones(1,N_active); %Receive apodization
% Rcv_apo = hann(N_active)';
[scat,start_time] = calc_scat_all(emit_aperture,receive_aperture,phantom_positions,phantom_amplitudes,1);
scat = [zeros(round(start_time*fs),N_active*N_active);scat];
%Apply transmit apodization
% rcv_apo_matrix = repmat(Rcv_apo,size(scat,1),1);
% rcv_apo_matrix = repmat(rcv_apo_matrix,1,N_active);
% trans_apo_matrix = zeros(size(rcv_apo_matrix));
% for i = 1:N_active
%     trans_apo_matrix(:,(i-1)*N_active+1:i*N_active) = Tx_apo(i);
% end
% scat_tx_apo = scat.*trans_apo_mat
% Reshape the synthetic aperture matrix
% N_transmit = N_elements;
% N_receive = N_elements;
% scat_reshape = zeros(size(scat,1),N_transmit,N_receive);
% for i = 1:N_receive
%     scat_reshape(:,:,i) = scat_final(:,i:N_receive:N_receive*N_transmit);
% end
% 
% rf_final = squeeze(sum(scat_reshape,2));
% figure;imagesc(sqrt(abs(rf_final)));

%%
time_samples = size(scat);
scat = reshape(scat,[time_samples(1),N_elements,N_elements]);
scat = sum(scat,2);
scat = squeeze(scat);

% scat = scat';
%% DAS beamforming
if sum(sum(sum(isnan(scat))))~= 0
    disp('NAN. Why?');
    u = u-1;
else
    disp('No problem.');
    disp(sum(sum(groundtruth)));
    Steer = deg2rad([-1,0,1]);  
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
    % prf = 500;                       % frame rate42
    PageEnd = size(RData,1);
    % Reshape the raw channel data to 3D (x,z,t)
    ReShpRData = zeros(PageEnd,chano);
    temp = zeros(Na,Nl,NoAngles);
    for g = 1:NoAngles
        Tx_focus = 0e-3; %m
        Tx_Delay = (arrayx'- Origin(1)) * sin(Steer(g))/c;
        Tx_Delay = Tx_Delay - min(Tx_Delay(:)); 
        dummy = double(squeeze(RData(1:PageEnd,:)));
        ReShpRData(:,:) = reshape(dummy,PageEnd,N_receive);

    %ReShpRData = scat for single plane wave imaging
    % Pixel-oriented beamforming (matrix version)

        clear delaysub delay delayidx;

        Tdelay = interp1(arrayx,Tx_Delay,x,'linear');
        Tdelay = repmat(reshape(Tdelay,[1,Nl,1]),[Na,1,chano]); % Transmit delay
%   delay = round(((Z+sqrt((Z-arrayZ).^2+(X-arrayX).^2))/c+Tdelay+lensCorrection/c*2+twpeak)*Fs);% calculate the delay in each channel for each pixel in units of samples
        delay = round(((Z * cos(Steer(g)) +sqrt((Z-arrayZ).^2+(X-arrayX).^2))/c+Tdelay+twpeak)*Fs);
        delay(delay<=0) = 1; % clear 0 delays
        delay(delay>PageEnd) = PageEnd; % clear delays beyond the sampling range
    %delaysub(:,:,:,n) = delay;
    %delayidx(:,:,:,n) = sub2ind([size(ReShpRData,1),chano],delay,repmat(reshape(1:chano,[1,1,chano]),[Na,Nl,1]));%find the indices of the data for each pixel in each channel
        delay = sub2ind([size(ReShpRData,1),chano],delay,repmat(reshape(1:chano,[1,1,chano]),[Na,Nl,1]));%find the indices of the data for each pixel in each channel
    % Beamforming and IQ downmixing
        dummy = squeeze(ReShpRData(:,:));
    %size(dummy) is 6547*128
        temp(:,:,g) = sum(dummy(delay),3);
    end
    %TEMP = sum(dummy(delayidx(:,:,:,n)).*ElementSens,3);
    %size(TEMP) is 1624*495, size(delayidx is 1624*495*128)
    %beamformedIQ(:,:,n,i) = TEMP;
    %env(:,:,n,i) = abs(hilbert(TEMP));
    X =[];
    Z =[];
    arrayZ =[];
    arrayX =[];
    Tdelay =[];
    %env = beamformedIQ;
    %data(:,:) = beamformedIQ;
    %figure;imagesc(x*1e3,z*1e3,20*log10(env./max(env(:))),[-60 0]);colormap(gray);
    %figure;imagesc(x*1e3,z*1e3,beamformedIQ);colormap(gray);
    beamformedIQ = zeros(Na,Nl);
    for v=1:NoAngles
    beamformedIQ = beamformedIQ+temp(:,:,v);
    end
    beamformedIQ = beamformedIQ/NoAngles;
    beamformedIQ = abs( hilbert(beamformedIQ));
    figure;imagesc(x*1e3,z*1e3,20*log10(beamformedIQ./max(beamformedIQ(:))),[-60 0]);colormap(gray);
    figure;imagesc(x*1e3,z*1e3,beamformedIQ);colormap(gray)
    
    beamformedIQ2 = temp(:,:,(NoAngles+1)/2);
    beamformedIQ2 = abs( hilbert(beamformedIQ2));
    figure;imagesc(x*1e3,z*1e3,20*log10(beamformedIQ2./max(beamformedIQ2(:))),[-60 0]);colormap(gray);
    figure;imagesc(x*1e3,z*1e3,beamformedIQ2);colormap(gray)
    %save(strcat('ahalf',int2str(u),'.mat'),'beamformedIQ');
    %save(strcat('ahalfgt',int2str(u),'.mat'),'groundtruth');
end
toc
end
