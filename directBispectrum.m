function Bspec = directBispectrum (X, fftLength, windowType, sampleSize, percentOverlap)
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


%% 2D Window Creation
if (exist('windowType', 'var') ~= 1) 
    windowType = 5; 
end
[m,n] = size(windowType);
window = windowType;

if (max(m,n) == 1)     % scalar: wind is size of Rao-Gabr window
    winsize = windowType;
    if (winsize < 0) 
        winsize = 5; 
    end        % the window length L
    winsize = winsize - rem(winsize,2) + 1;  % make it odd
    if (winsize > 1)
        mwind   = fix (fftLength/winsize);            % the scale parameter M
        lby2    = (winsize - 1)/2;

        theta  = -lby2:lby2;
        opwind = ones(winsize,1) * (theta .^2);       % w(m,n)=m^2
        opwind = opwind + opwind' + theta' * theta;   % m^2 + n^2 + mn
        opwind = 1 - (2*mwind/fftLength)^2 * opwind;       %
        hex    = ones(winsize,1) * theta;             % m
        hex    = abs(hex) + abs(hex') + abs(hex+hex');
        hex    = (hex < winsize);
        opwind = opwind .* hex;
        opwind = opwind * (4 * mwind^2) / (7 * pi^2) ;
    else
        opwind = 1;
    end

elseif (min(m,n) == 1)  % 1-D window passed: convert to 2-D
    window = window(:);
    if (any(imag(window) ~= 0))
        disp('1-D window has imaginary components: window ignored')
        window = 1;
    end
    if (any(window < 0))
        disp('1-D window has negative components: window ignored')
        window = 1;
    end
    lwind  = length(window);
    windf  = [window(lwind:-1:2); window];  % the full symmetric 1-D
    window = [window; zeros(lwind-1,1)];
    opwind = (windf * windf').* hankel(flipud(window), window); % w(m)w(n)w(m+n)
    winsize = length(window);

else                    % 2-D window passed: use directly
    winsize = m;
    if (m ~= n)
        disp('2-D window is not square: window ignored')
        window = 1;
        winsize = m;
    end
    if (rem(m,2) == 0)
        disp('2-D window does not have odd length: window ignored')
        window = 1;
        winsize = m;
    end
    opwind  = window;
end

%% Accumulate Triple Products
Bspec    = zeros(fftLength,fftLength);

mask = hankel(1:fftLength,[fftLength,1:fftLength-1] );   % the hankel mask (faster)
locseg = (1:sampleSize)';
for krec = 1:numCols
    xseg   = X(locseg);
    Xf     = fft(xseg-mean(xseg), fftLength)/sampleSize;
    CXf    = conj(Xf);
    Bspec  = Bspec + (Xf * Xf.') .*reshape(CXf(mask), fftLength, fftLength);
    locseg = locseg + nAdvance;
end

Bspec = fftshift(Bspec)/(numCols);

%% Frequency-Domain Smoothing
if (winsize > 1)
    lby2 = (winsize-1)/2;
    Bspec = conv2(Bspec,opwind);
    Bspec = Bspec(lby2+1:lby2+fftLength,lby2+1:lby2+fftLength);
end


%% Contour Plot of Magnitude Bispectum
if (rem(fftLength,2) == 0)
    waxis = (-fftLength/2:(fftLength/2-1))'/fftLength;
else
    waxis = (-(fftLength-1)/2:(fftLength-1)/2)'/fftLength;
end

contour(waxis,waxis,abs(Bspec),4),grid on 
title('Bispectrum estimated via the direct (FFT) method')
xlabel('f1'), ylabel('f2')
end