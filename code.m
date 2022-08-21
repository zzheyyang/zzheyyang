%%-----ECG------
clear, clc, close all
[tm,data1]=rdsamp('01',1,13000); %read data
figure();
plot(tm);
xlabel('Samples');ylabel('Amplitude'); axis([12000 13000 -1 1]);grid on; box on;

x=tm;
fs=128;

ecg=tm;
if ~isvector(ecg)
  error('ecg must be a row or column vector');
end
ecg = ecg(:); % vectorize


Fs = 128;
flag = 1;
%Ecg is the input signal, Fs is the sampling frequency, plot is whether to plot, plot is 1 and plot is 0


% The signal and threshold after filtering to find the derivative square mean
Signal_mean = [];
Signal_mean_index = [];
Thresh_Signal_mean = 0;
Thresh_Noise_mean = 0;

% Noise after filtering to find the derivative square mean
Noise_mean = [];
Noise_mean_index = [];

% Signal and threshold after band pass filtering
Signal_filter = [];
Signal_filter_index = [];
Thresh_Signal_filter = 0;
Thresh_Noise_filter = 0;

% Current signal, current noise level
Signal_temp_mean = 0;
Noise_temp_mean = 0;
Signal_temp_filter = 0;
Noise_temp_filter = 0;

% Transient storage of signal and noise after filtering to find the derivative square mean
Signal_mean_buf = [];
Noise_mean_buf = [];


% Transient storage of signal and noise after band pass filtering
Signal_filter_buf = [];
Noise_filter_buf = [];

Thresh_Signal_mean_buf = [];
Thresh_Signal_filter_buf = [];

% 8 regular heart rate intervals
regular_RR_interval = 0;
% 8 heart rate intervals
mean_interval = 0;
% Heart rate interval, regular_RR_interval is regular_RR_interval, otherwise it is mean_interval
RR_interval =  0;

%When <0.36Fs, determine whether it is T wave. If it is indeed R wave, it is 0; otherwise, it is 1
R_Flag = 0;

%% Filter 1, low pass then high pass
if Fs == 128
    
    %low-pass，H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
    
    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a = [1 -2 1];
    
    h_l = filter(b,a, [1,zeros(1,12)]);
    
    ecg_l = conv(ecg, h_l);
    ecg_l = ecg_l / max(abs(ecg_l));
    
    
    %high pass H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
    
    b = [-1 zeros(1,15) 32 -32 zeros(1,14) 1];
    a = [1 -1];
    
    h_h = filter(b,a, [1, zeros(1,32)]);
    ecg_filter = conv(ecg_l, h_h);
    ecg_filter = ecg_filter / max(abs(ecg_filter));

%% Filter 2, Direct Buttord filter, 5-15Hz digital bandpass filter
else
    wp = [5 15] * 2 / Fs;
    
    [B,A] = butter(3, wp);%3-tap
    
    ecg_filter = filtfilt(B,A, ecg);
    ecg_filter = ecg_filter / max(abs(ecg_filter));
    
end

%% Derivative H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
    
    h_d = [-1 -2 0 2 1]/8;
    ecg_deri = conv(h_d, ecg_filter);
    ecg_deri = ecg_deri / max(abs(ecg_deri));

%% Square
    
    ecg_square = ecg_deri .^ 2;
    
%% Window averaging
    Window_width = round(0.15 * Fs);
    
    ecg_mean = conv(ecg_square, ones(1, Window_width) / Window_width);

%% Plot
if flag
    figure;
    if Fs == 128
    ax(1)=subplot(321);plot(ecg);axis tight;title('Original ECG signal');
    axis([6000 7000 -1 1]);grid on; box on;xlabel('Samples');ylabel('Amplitude')
    ax(2)=subplot(322);plot(ecg_l);axis tight;title('Low pass filtered signal');
    axis([6000 7000 -1 1]);grid on; box on;xlabel('Samples');ylabel('Amplitude')
    else
        ax(1)=subplot(3,2,[1,2]);plot(ecg);axis tight;title('Original ECG signal');
        axis([6000 7000 -1 1]);grid on; box on;xlabel('Samples');ylabel('Amplitude')
    end
    
    ax(3) = subplot(323);plot(ecg_filter);axis tight;title('Band pass filtered signal');
    axis([6000 7000 -1 1]);grid on; box on;xlabel('Samples');ylabel('Amplitude')
    ax(4) = subplot(324);plot(ecg_deri);axis tight;title('Differential operation after signal');
    axis([6000 7000 -1 1]);grid on; box on;xlabel('Samples');ylabel('Amplitude')
    ax(5) = subplot(325);plot(ecg_square);axis tight;title('Squared signal');
    axis([6000 7000 -1 1]);grid on; box on;xlabel('Samples');ylabel('Amplitude')
    ax(6) = subplot(326);plot(ecg_mean);axis tight;title('The sliding window detects the signal');
    axis([3000 7000 -1 1]);grid on; box on;xlabel('Samples');ylabel('Amplitude')
    
end


%% Find local peaks by findpeaks

[peaks, locs] = findpeaks(ecg_mean, 'MINPEAKDISTANCE', round(0.2 * Fs));


%% Initialize the threshold using the data from the first 2s


% Initialize signal thresholds and noise thresholds using only square and mean the first two seconds after the signal
% The signal threshold is 0.33 of the maximum value of the first two seconds, and the noise threshold is 0.5 of the average value of the first two seconds
Thresh_Signal_mean = max(ecg_mean(1:2*Fs)) / 3;
Thresh_Noise_mean = mean(ecg_mean(1:2*Fs)) / 2;

Signal_temp_mean = Thresh_Signal_mean;
Noise_temp_mean = Thresh_Noise_mean;

% bandpass filtering signal
%
Thresh_Signal_filter = max(ecg_filter(1:2*Fs)) / 3;
Thresh_Noise_filter = mean(ecg_filter(1:2*Fs)) / 2;

Signal_temp_filter = Thresh_Signal_filter;
Noise_temp_filter = Thresh_Noise_filter;


%% Each peak is treated
for i = 1 : length(peaks)  
    
    %% Find the position of the signal peak after filtering, index is relative value
    if locs(i) - round(0.15 * Fs) >= 1 && locs(i) <= length(ecg_filter)
        [peak_filter, index_filter] = max( ecg_filter(locs(i) - round(0.15 * Fs) : locs(i)));
        index_filter = locs(i) - round(0.15 * Fs) + index_filter - 1;
    else
        if i == 1 % i = 1 is when index is absolute value
            [peak_filter, index_filter] = max( ecg_filter(1 : locs(i)));
        else %index is relative value
            [peak_filter, index_filter] = max( ecg_filter(locs(i) - round(0.15 * Fs) : end));
            index_filter = locs(i) - round(0.15 * Fs) + index_filter - 1 ;
        end
    end
    %% If the interval is greater than 1.66, the missing R wave should be processed first, so as to maintain the order
    
    % It takes 9 R waves to have 8 intervals to determine if it is greater than 1.66
    if length(Signal_mean) >= 9
        RR_diff = diff(Signal_mean_index(end - 8:end));
        mean_interval = mean(RR_diff);
        
        temp_interval = Signal_mean_index(end) - Signal_mean_index(end - 1);
        
        if(temp_interval >= 0.92 * mean_interval || temp_interval <= 1.16 * mean_interval)
            Thresh_Signal_mean = Thresh_Signal_mean / 2;
            Thresh_Signal_filter = Thresh_Signal_filter / 2;
        else
            regular_RR_interval = mean_interval;
        end
    end
    
    if regular_RR_interval
        RR_interval = regular_RR_interval;
    elseif mean_interval
        RR_interval = mean_interval;
    else
        RR_interval = 0;
    end
    
    % If the heart rate interval has value
    if RR_interval
        % More than 1.66 RR, indicating that there is a missed detection, refractory period is 0.2s
        if locs(i) - Signal_mean_index(end) >= 1.66 * RR_interval
            [pk_temp, pk_temp_ind] = max( ecg_mean(Signal_mean_index(end) + round(0.2 * Fs) : locs(i) - round(0.2 * Fs)));
            
            pk_temp_ind = Signal_mean_index(end) + round(0.2 * Fs) + pk_temp_ind - 1;
            
            % The mean signal is greater than the noise peak and is added to buffer
            if pk_temp >= Thresh_Noise_mean
                Signal_mean = [Signal_mean pk_temp];
                Signal_mean_index = [Signal_mean_index pk_temp_ind];
                
                %Find the location of the filter signal
                if pk_temp_ind <= length(ecg_filter)
                    [pk_filter_temp, pk_filter_temp_ind] = max(ecg_filter(pk_temp_ind - round(0.15 * Fs) :  pk_temp_ind));
                else
                    [pk_filter_temp, pk_filter_temp_ind] = max(ecg_filter(pk_temp_ind - round(0.15 * Fs) :  length(ecg_filter)));
                end
                
                pk_filter_temp_ind = ecg_filter(pk_temp_ind - round(0.15 * Fs) + pk_filter_temp_ind - 1);
                
                % The filtered signal is also greater than the noise threshold and added to the BUF, while the SIG uses a new update strategy
                if pk_filter_temp >= Thresh_Noise_filter
                    Signal_filter = [Signal_filter pk_filter_temp];
                    Signal_filter_index = [Signal_filter_index pk_filter_temp_ind];
                    
                    Signal_temp_filter = 0.25 * pk_filter_temp + 0.75 * Signal_temp_filter;
                end
                    
                Signal_temp_mean = 0.25 * pk_temp + 0.75 * Signal_temp_mean;

            end
        end
    end
    
    %% The current wave peak exceeds the signal threshold
    if peaks(i) >= Thresh_Signal_mean
        
        % have three waves to have two intervals
        if(length(Signal_mean) >= 3)
            
            %If the interval is less than 0.36fs, determine whether it is T wave or noise
            if locs(i) - Signal_mean_index(end) < round(0.36 * Fs)
                Slope_1 = mean(diff(ecg_mean(locs(i) - round(0.075 * Fs) : locs(i))));
                Slope_2 = mean(diff(ecg_mean(Signal_mean_index(end) - round(0.075 * Fs) : Signal_mean_index(end))));
                
                %The slope is too small and the current crest is identified as noise
                if Slope_1 <= Slope_2 / 2
                    
                    Noise_mean = [Noise_mean peaks(i)];
                    Noise_mean_index = [Noise_mean_index locs(i)];
                    Noise_temp_mean = 0.25 * peaks(i) + 0.75 * Noise_temp_mean;
                    Noise_temp_filter = 0.25 * peak_filter + 0.75 * Noise_temp_filter;
                    R_Flag = 1;%is T wave
                
                else %The slope is high enough. It's an R wave
                    R_Flag = 0;%T wave
                end
            end
        end
        
        if R_Flag == 0

            Signal_mean = [Signal_mean peaks(i)];
            Signal_mean_index = [Signal_mean_index locs(i)];
            Signal_temp_mean = 0.125 * peaks(i) + 0.875 * Signal_temp_mean;

            if peak_filter >= Thresh_Signal_filter

                Signal_filter = [Signal_filter peak_filter];
                Signal_filter_index = [Signal_filter_index index_filter];
                Signal_temp_filter = 0.125 * peak_filter + 0.875 * Signal_temp_filter;

            end
        end
        
                    
    %% The current wave peak is greater than the noise threshold, but less than the signal threshold, and the wave peak is the noise peak
    elseif peaks(i) < Thresh_Signal_mean && peaks(i) >= Thresh_Noise_mean
        
        Noise_temp_mean = 0.125 * peaks(i) + 0.875 * Noise_temp_mean;
        Noise_temp_filter = 0.125 * peak_filter + 0.875 * Noise_temp_filter;
        
    %% less than,noise peak
    else
        Noise_temp_mean = 0.125 * peaks(i) + 0.875 * Noise_temp_mean;
        Noise_temp_filter = 0.125 * peak_filter + 0.875 * Noise_temp_filter;
        
        Noise_mean = [Noise_mean peaks(i)];
        Noise_mean_index = [Noise_mean_index locs(i)];
        
    end
    
    %% Update threshold parameters
    Thresh_Signal_mean = Noise_temp_mean + 0.25 * (Signal_temp_mean - Noise_temp_mean);
    Thresh_Noise_mean = Thresh_Signal_mean / 2;
    
    Thresh_Signal_filter = Noise_temp_filter + 0.25 * (Signal_temp_filter - Noise_temp_filter);
    Thresh_Noise_filter = Thresh_Signal_filter / 2;
    
    %% The current value of the signal is stored in BUF
    Signal_mean_buf = [Signal_mean_buf Signal_temp_mean];
    Noise_mean_buf = [Noise_mean_buf Noise_temp_mean];
    

    % Transient storage of signal and noise after band pass filtering
    Signal_filter_buf = [Signal_filter_buf Signal_temp_filter];
    Noise_filter_buf = [Noise_filter_buf Noise_temp_filter];
    
    Thresh_Signal_mean_buf = [Thresh_Signal_mean_buf Thresh_Signal_mean];
    Thresh_Signal_filter_buf = [Thresh_Signal_filter_buf Thresh_Signal_filter];
    
    
    % reset Flag
    R_Flag = 0;
    
end


%% plot
if flag
    figure;
    plot(ecg_filter);axis tight;title('R wave location results');axis([6000 7000 -1 1])%axis([12000 13000 -1 1]);
    grid on; box on;xlabel('Samples');ylabel('Amplitude')
    hold on
    scatter(Signal_filter_index, Signal_filter);
end

%distance=Signal_filter_index(60)-Signal_filter_index(10);%long term AF,fs=128hz
distance=max(Signal_filter_index)-Signal_filter_index(5);%fs=300hz
[m,p]=max(Signal_filter_index);
t=distance/fs;
heart_rate=(p-5-1)./t*60;
fprintf('HR(bmp)=%f\n',heart_rate);
