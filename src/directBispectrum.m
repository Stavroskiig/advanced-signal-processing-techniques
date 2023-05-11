function bispectrum = directBispectrum (X, fftLength, windowType, sampleSize, percentOverlap)
%% Parameter Validation

% Get the number of rows and columns of the input data matrix X
[numRows, numCols] = size(X);

% If X has only one row, reshape it into a column vector
if (numRows == 1) 
    X = X(:);  
    numRows = numCols; 
    numCols = 1; 
end

% Check if the variable 'fftLength' is not defined, set it to the default value 128
if (exist('fftLength','var') ~= 1)            
    fftLength = 128; 
end

% Check if the variable 'percentOverlap' is not defined, set it to the default value 50
if (exist('percentOverlap','var') ~= 1)      
    percentOverlap = 50;  
end

% Ensure that percentOverlap is within the range of 0 to 99
percentOverlap = min(99, max(percentOverlap, 0));

% If there are multiple columns in X, set percentOverlap to 0
if (numCols > 1)                  
    percentOverlap =  0;  
end

% Check if the variable 'sampleSize' is not defined, set it to 0
if (exist('sampleSize','var') ~= 1)          
    sampleSize = 0;   
end

% If there are multiple columns in X, set sampleSize to the number of rows
if (numCols > 1)                    
    sampleSize = numRows;  
end

% If X has only one column and sampleSize is not positive, compute sampleSize based on percentOverlap
if (numCols == 1 && sampleSize <= 0)
sampleSize = fix(numRows/ (8 - 7 * percentOverlap/100));
end

% If fftLength is smaller than sampleSize, set it to the smallest power of 2 that is greater than or equal to sampleSize
if (fftLength  < sampleSize)   
    fftLength = 2^nextpow2(sampleSize); 
end

% Compute the number of columns based on numRows, numCols, percentOverlap, and nadvance
percentOverlap = fix(sampleSize * percentOverlap / 100);
nAdvance = sampleSize - percentOverlap;
numCols = fix( (numRows*numCols - percentOverlap) / nAdvance);


%% Window Creation

% Check if the windowType variable exists, otherwise assign a default value
if (exist('windowType', 'var') ~= 1) 
    windowType = 5; 
end

% Get the size of the windowType variable
[winRows, winCols] = size(windowType);

% Assign the window variable based on the windowType value
window = windowType;

% Check if the windowType is a scalar (indicating the size of Rao-Gabr window)
if (max(winRows, winCols) == 1)
    % Set the window size based on the windowType value
    windowSize = windowType;

    % Check if winsize is less than 0, assign a default value
    if (windowSize < 0) 
        windowSize = 5; 
    end

    % Make winsize odd
    windowSize = windowSize - rem(windowSize,2) + 1;

    % Perform window creation if winsize is greater than 1
    if (windowSize > 1)
        scaleParam   = fix (fftLength / windowSize); % Scale parameter M
        halfWindowSize    = (windowSize - 1) / 2;

        % Create the window coefficients based on mathematical operations
        theta  = -halfWindowSize:halfWindowSize;
        opwind = ones(windowSize,1) * (theta .^2); % w(m,n)=m^2
        opwind = opwind + opwind' + theta' * theta; % m^2 + n^2 + mn
        opwind = 1 - (2*scaleParam/fftLength)^2 * opwind; 
        hex    = ones(windowSize,1) * theta; 
        hex    = abs(hex) + abs(hex') + abs(hex+hex');
        hex    = (hex < windowSize);
        opwind = opwind .* hex;
        opwind = opwind * (4 * scaleParam^2) / (7 * pi^2) ;
    else
        % Window size is 1, assign 1 as the window coefficient
        opwind = 1;
    end

% Check if the windowType is a 1D window passed as input
elseif (min(winRows, winCols) == 1)
    % Convert the 1D window to a 2D window
    window = window(:);

    % Check if the 1D window has imaginary or negative components, if so, assign 1 as the window
    if (any(imag(window) ~= 0))
        disp('1D window has imaginary components: window ignored')
        window = 1;
    end
    if (any(window < 0))
        disp('1D window has negative components: window ignored')
        window = 1;
    end
    windowLength  = length(window);
    fullSymmetricWindow  = [window(windowLength:-1:2); window];  % the full symmetric 1-D
    window = [window; zeros(windowLength-1,1)];
    opwind = (fullSymmetricWindow * fullSymmetricWindow').* hankel(flipud(window), window); % w(m)w(n)w(m+n)
    windowSize = length(window);

% Check if the windowType is a 2D window passed as input
else
    windowSize = winRows;

    % Check if the 2D window is square and has odd length, otherwise assign 1 as the window
    if (winRows ~= winCols)
        disp('2D window is not square: window ignored')
        window = 1;
        windowSize = numRows;
    end
    if (rem(winRows,2) == 0)
        disp('2D window does not have odd length: window ignored')
        window = 1;
        windowSize = winRows;
    end

    opwind  = window;
end

%% Accumulate Triple Products

% Initialize the bispectrum matrix
bispectrum = zeros(fftLength, fftLength);

% Create a hankel mask for faster processing
hankelMask = hankel(1:fftLength, [fftLength, 1:fftLength-1]);

% Initialize the segment location
segmentLocation = (1:sampleSize)';

% Iterate over each column of the input data
for columnIdx = 1:numCols
    % Extract the segment of the input data
    segment = X(segmentLocation);

    % Compute the Fourier transform of the segment (normalized)
    segmentFourier = fft(segment - mean(segment), fftLength) / sampleSize;

    % Compute the complex conjugate of the Fourier transform
    complexConjugate = conj(segmentFourier);

    % Compute the outer product of the Fourier transform and its conjugate
    % and accumulate it in the bispectrum matrix
    bispectrum = bispectrum + (segmentFourier * segmentFourier.') .* reshape(complexConjugate(hankelMask), fftLength, fftLength);

    % Update the segment location by advancing it
    segmentLocation = segmentLocation + nAdvance;
end

% Normalize the bispectrum by dividing it by the number of columns
bispectrum = fftshift(bispectrum) / numCols;

%% Frequency-Domain Smoothing

% Check if windowSize is greater than 1 for smoothing
if (windowSize > 1)
    halfWindowSize = (windowSize-1) / 2;

    % Apply convolution with the window in the frequency domain
    bispectrum = conv2(bispectrum, opwind);
    
    % Extract the central part of the smoothed bispectrum after convolution
    bispectrum = bispectrum(halfWindowSize+1:halfWindowSize+fftLength,halfWindowSize+1:halfWindowSize+fftLength);
end

%% Contour Plot of Magnitude Bispectum

% Compute the frequency axis values
if (rem(fftLength,2) == 0)
    % For even fftLength, use the appropriate frequency axis values
    waxis = (-fftLength/2:(fftLength/2-1))'/fftLength;
else
    % For odd fftLength, use the appropriate frequency axis values
    waxis = (-(fftLength-1)/2:(fftLength-1)/2)'/fftLength;
end

% Plot the magnitude of the bispectrum using contour plot
contour(waxis,waxis,abs(bispectrum),4),grid on 

% Set the title and labels for the plot
%title('Bispectrum Estimation using the Direct Method')
xlabel('f1'), ylabel('f2')
end