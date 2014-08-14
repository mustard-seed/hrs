function jeppler_hsrsar_backproj
%physical constants
c = 299792458.0;

%define image formation parameters
dx = 1; dy=2; dz=1;
x = -10:dx:100;
y = 0:dy:0;
z = 0:dz:100;
doppler_fraction = 0.9; %fraction of full Doppler spectrum to be used for each target - small fraction means lower resolution but faster processing

%define system parameters
fc = 9e8;
fc_band = 30.0e6;
del_fc = 200e3;
fd_samp = 1000.0; %baseband sample rate (Hz)
t_range = [0, 10.0]; %data acquisition time interval

%define geometry
p_t = [0, 0, 30]; %transmitter location
v = [0.0, 100.0, 0.0]; %constant velocity of receiver on train
p_r_0 = [50, -500, 4.5]; %initial position of receiver on train
p_s_1 = [0, 0, 0.0]; %scatterer 1

%derived values
t = t_range(1): 1/fd_samp:t_range(2);
fc_vec = fc-fc_band/2:del_fc:fc+fc_band/2;
fc_band = fc_vec(end) - fc_vec(1);
img = zeros(length(x), length(y), length(z));
img = complex(img, img);
alpha=0.0; %random phase offset

%simulate data acquisition
signal = get_signal(p_t, p_s_1,fc_vec,alpha,v, p_r_0, t); %baseband signal per tone

%range compression
signal_rc = fft(signal,[],1); %range comressed signal
signal_rc_ex = [signal_rc' (signal_rc(1,:))']'; %concatenate zero range row to the end to enable wrapped interpolation

%JDE_TODO correct range compressed data for free space path loss
FSPL = 20 * log10(rTotal)+ 20 * log10(fc(i_fc))-147.55;


%set-up interpolation of range compressed data
[n_rng, n_az] = size(signal_rc);
del_rng = c/fc_band; %m
ranges = (0:n_rng-1)*del_rng;
[X, Y] = ndgrid([ranges, ranges(end)+del_rng], t);
GI = griddedInterpolant(X,Y, signal_rc_ex); %set-up interpolation object

%back-projection over image grid
range_interval = n_rng*del_rng;
lambda = c/fc_vec(1);
del_y_factor = (doppler_fraction/sqrt(1-doppler_fraction^2));
for ii = 1:length(x)
    ii
    for jj = 1:length(y)
        for kk = 1:length(z)
            p_s = [x(ii), y(jj), z(kk)];
            del_xz = sqrt((p_r_0(1) - p_s(1))^2 + (p_r_0(3) - p_s(3))^2); %distance from current scatterer to Rx (assumes y only train velocity)
            del_y = del_y_factor*del_xz;
            t_range_segment = [max((p_s(2) - del_y - p_r_0(2))/v(2), t_range(1)), min((p_s(2) + del_y - p_r_0(2))/v(2), t_range(2))]; %again, train velocity assumed here
            t_segment = t_range_segment(1): 1/fd_samp:t_range_segment(2);
            rTotal_ = r_total(p_t, p_s,v, p_r_0, t_segment);
            rTotal = mod(rTotal_, range_interval);
            img(ii,jj,kk) = sum(GI(rTotal, t_segment).*exp(complex(0,-2*pi*(rTotal)./lambda)));
        end
    end
end

%render results for xy, xz, yz planes centered on scatterer
close all;
p_s_1_idx = round((p_s_1 - [x(1), y(1), z(1)])./[dx, dy, dz])+1;
%figure('units','normalized','outerposition',[0 0 1 1]);
figure;
subplot(2,2,1);
imagesc(x,y,abs(squeeze(img(:,:,p_s_1_idx(3)))'));
set(gca,'YDir','normal');
colormap(gray);

subplot(2,2,2);
imagesc(x,z,abs(squeeze(img(:,p_s_1_idx(2),:))'));
set(gca,'YDir','normal');
colormap(gray);

subplot(2,2,3);
imagesc(y,z,abs(squeeze(img(p_s_1_idx(1),:,:))'));
set(gca,'YDir','normal');
colormap(gray);

figure;
imagesc(x,z,abs(squeeze(img(:,p_s_1_idx(2),:))'));
set(gca,'YDir','normal');
colormap(gray);

img_norm = img./max(max(max(abs(img))));
%figure('units','normalized','outerposition',[0 0 1 1]);
%subplot(2,1,1);
%plot(20*log10(abs(img_norm(51,:))));
%subplot(2,1,2);
%plot(20*log10(abs(img_norm(:,51))));
end

%--------------------------------------------------------------------------
function rTotal = r_total(p_t_, p_s_,v, p_r_0, t)
    %Tx, scatterer and Rx positions
    nt = length(t);
    p_t = p_t_'*ones(1,nt);
    p_s = p_s_'*ones(1,nt);
    p_r = p_r_0'*ones(1,nt) + v'*t;

    %two leg vector
    Rts = p_s - p_t;
    Rsr = p_r - p_s;

    %two leg ranges
    Rts_norm = sqrt(Rts(1,:).^2+Rts(2,:).^2+Rts(3,:).^2);
    Rsr_norm = sqrt(Rsr(1,:).^2+Rsr(2,:).^2+Rsr(3,:).^2);
    
    %total range
    rTotal = Rts_norm + Rsr_norm;
end
%--------------------------------------------------------------------------

function FSPL = FSPL(range, f)
    FSPL_log = 20 * log10(range)+ 20 * log10(f)-147.55; 
    FSPL = 10.^(-FSPL_log/20);
end

%--------------------------------------------------------------------------
function signal = get_signal(p_t, p_s, fc, alpha, v, p_r_0, t)
c = 299792458.0;
rTotal = r_total(p_t, p_s, v, p_r_0, t);
signal = zeros(length(fc), length(rTotal));
for i_fc = 1:length(fc)
    lambda = c/fc(i_fc);
    
    % Free space path loss
    %FSPL = 20 * log10(rTotal)+ 20 * log10(fc(i_fc))-147.55; 
    %FSPL_linear = 10.^(-FSPL/20);
    
    %JDE_TODO undo this!!!
    FSPL = FSPL(rTotal,fc(i_fc));
    %FSPL = ones(1, length(rTotal));
    
    % Receive Power
    signal(i_fc, :) = FSPL.*exp(complex(0,alpha + 2*pi*(rTotal)./lambda));
end
end