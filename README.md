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
4. Range-Doppler Map(2D-FFT) : 1D-fft was performed on the single half column elements of the mix to show that range from 1d fft can be calculated if it is performed in the direction of changing rows. After that 2D-FFT is performed on the whole matrix which is formed from the mix signals/vectors matrix and is knowns as Range-Doppler Map(RDM).
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
5. CFAR Detection: A 2D CFAR is calculated based on the RDM in step 4. CFAR is calculated in such a way that edges of the matrix are not used to determine to CUT. Here Tr means training cells along range direction, Gr means guard cells along the range direction, Td means training cells along the doppler direction and Gd means guard cells along the doppler direction. 
   ```
   for i = Tr+Gr+1:(Nr/2)-(Gr+Tr)  %Nr/2 because of taking one side of the spectrum see line 127 and 128
    for j = Td+Gd+1:Nd-(Gd+Td)
        noise_level = zeros(1,1);
        
        for p = i-(Tr+Gr) : i+(Tr+Gr)
            for q = j-(Td+Gd) : j+(Td+Gd)
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        
        threshold = pow2db(noise_level / (2 * (Tr+Gr+1) * 2 * (Td+Gd+1) - (Gr * Gd)-1)); %pow2db on normalized noise level 
        threshold = threshold + Offset;
        CUT = RDM(i, j);
        
        if (CUT < threshold)
            RDM(i, j) = 0;
        else
            RDM(i, j) = 1;
        end
        
      end
   end
   ```    
## 2. Selection of Training, Guard cells and offset

```
    %Select the number of Training Cells in both the dimensions.
    Tr = 12; %Training Cells in range dimension meaning along changing cols dir.
    Td = 10; %Training Cells in doppler dimension meaning along changing rows dir.

    %Select the number of Guard Cells in both dimensions around the Cell under 
    %test (CUT) for accurate estimation
    Gr = 4;
    Gd = 4;

    % offset the threshold by SNR value in dB
    Offset = 1.2; 
```

## 3. Steps taken to suppress the non-thresholded cells at the edges
As CUT can not be located at the edges therefore cells on the edges can't take part in the thresholding process. Also to detect the signal in CUT, all values below the average noise floor(final threshold) of matrix elements which were also less than the threshold were set to 0 otherwise 1. This means we can set all values in matrix to zero where values are not zero or not 1. So in the end only those values remain in the RDM which were set to 1 by the 2D-CFAR algorithm.
```
    RDM(RDM ~= 0 & RDM ~= 1) = 0;

```