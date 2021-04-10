# SFND_Radar_Target_Generation_and_Detection
![](CFAR_2D_Process.png "Figure 1: Block diagram of the entire final project to compute 2D-CFAR")

## 1. Implementation steps for the 2D CFAR process
From Figure 1. the following steps can be concluded and implemented the final project:
1. FMCW Configuration: From the given FMCW radar parameters a specific configuration is generated using bandwidth, chirp time and slope formulas:
 ```
    Fc= 77e9; %carrier freq
    MAX_RANGE = 200;
    MAX_VELOCITY = 70; %m/s
    RANGE_RESOLUTION = 1; %1meter
    C = 3e8;

    B_SWEEP = C/2 * RANGE_RESOLUTION;
    T_CHIRP = 5.5 * 2 * MAX_RANGE / C;
    SLOPE = B_SWEEP / T_CHIRP; 

 ```
 2. Moving Target Generation: Based on the initial position and velocity a simulation of based on the constant velocity model is generated. 
 ```
    INITIAL_RANGE = 110; %meters initial range
    INITIAL_VELOCITY = -20; %m/s initial velocity/speed of the target
    for i=1:length(time_vector)         
    
        %For each time stamp update the Range of the Target for constant velocity. 
        range_covered(i) = INITIAL_RANGE + INITIAL_VELOCITY * time_vector(i);
        time_delay(i) = 2 * range_covered(i) / C;
        ....
    end    
 ```
3. Signal Propagation and Processing Reflected Signal: Tx and Rx signals of a moving target are modeled based on the frequency modulation of the radar parameter such as slope and time delay. Here it can be seen that Rx signal is just time delayed version of the Tx signal. Beat frequencies or the mix signals are calculated by mutliplying Tx and Rx vectors in such a way that multiplication produce differences in frequency components where larger difference means the object is farther from the ego vehicle and vice versa:
   ```
   for i=1:length(time_vector)         
    
        %For each time stamp update the Range of the Target for constant velocity. 
        range_covered(i) = INITIAL_RANGE + INITIAL_VELOCITY * time_vector(i);
        time_delay(i) = 2 * range_covered(i) / C;
    
        %For each time sample we need update the transmitted and
        %received signal.
        % FM of CW takes place in Tx and Rx
        Tx(i) = cos(2 * pi * (Fc * time_vector(i) + SLOPE * time_vector(i)^2/2));
        Rx(i) = cos(2 * pi * (Fc * (time_vector(i)-time_delay(i)) + SLOPE * (time_vector(i)-time_delay(i))^2/2));
        
        %Now by mixing the Transmit and Receive generate the beat signal
        %This is done by element wise matrix multiplication of Transmit and
        %Receiver Signal
        Mix(i) = Tx(i) * Rx(i);
    
    end
   ```  
4. Range-Doppler Map(2D-FFT) : 1D-fft was performed on the single half column elements of the mix to show that range from 1d fft can be calculated if it is performed in the direction of chaning rows. After that 2D-FFT is performed on the whole matrix which is formed from the mix singls and is knowns as Range-Doppler Map(RDM).
   ```
    Mix_reshape = reshape(Mix, Nr, Nd);

    %run the FFT on the beat signal along the range bins dimension (Nr) and
    fft_range = fft(Mix_reshape, Nr);

    % Take the absolute value of FFT output
    fft_range_abs = abs(fft_range);

    % Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
    % Hence we throw out half of the samples.
    fft_range_normalized = fft_range_abs ./ max(fft_range_abs) ; %elemewise normalization 
    %fft_range_half_spectrum = fft_range_normalized(1:Nr/2, 1:2);
    fft_range_half_spectrum = fft_range_normalized(1:Nr/2); %beat signal in a single half col. 

    %plotting the range
    figure ('Name', 'Range from First FFT(Range that is extracted from single col.)')
    subplot(2,1,1)

    % plot FFT output 
    plot(fft_range_half_spectrum); 
    axis ([0 200 0 1]);


    Mix = reshape(Mix, [Nr,Nd]);
    % 2D FFT using the FFT size for both dimensions.
    sig_fft2 = fft2(Mix, Nr, Nd);

    % Taking just one side of signal from Range dimension.
    sig_fft2 = sig_fft2(1:Nr/2, 1:Nd);
    sig_fft2 = fftshift (sig_fft2);
    RDM = abs(sig_fft2);
    RDM = 10 * log10(RDM) ;

    %use the surf function to plot the output of 2DFFT and to show axis in both
    %dimensions
    doppler_axis = linspace(-100, 100, Nd);
    range_axis = linspace(-200,200,Nr/2) * ((Nr/2)/400);
    figure('Name', 'Range-Doppler Map');
    surf(doppler_axis, range_axis, RDM);
    xlabel('Velocity');
    ylabel('Range');
    zlabel('RDM Amplitude')
   ```
## 2. Selection of Training, Guard cells and offset


## 3. Steps taken to suppress the non-thresholded cells at the edges