function bispectrum = indirectBispectrum (X, maxLag, sampleSize,percentOverlap, biasCorrection, fftLength, windowType)
%% Parameter Validation

% Get the number of rows and columns of the input data matrix X
[numRows, numCols] = size(X); 

% If X has only one row, reshape it into a column vector
if (numRows == 1) 
    X = X(:);   
    numRows = numCols; 
    numCols = 1;      
end 

% Check if the variable 'percentOverlap' is not defined, set it to the default value 50
if (exist('percentOverlap','var') ~= 1)   
    percentOverlap = 0;          
end

% Ensure that percentOverlap is within the range of 0 to 99
percentOverlap = min(99, max(percentOverlap, 0)); 

% If there are multiple columns in X, set percentOverlap to 0
if (numCols > 1)               
    percentOverlap = 0;          
end 

% Check if the variable 'sampleSize' is not defined, set it to numRows
if (exist('sampleSize', 'var') ~= 1)     
    sampleSize = numRows;         
end

% Ensure that sampleSize is within the valid range
if (sampleSize > numRows || sampleSize <= 0) 
    sampleSize = numRows;         
end

% Check if the variable 'biasCorrection' is not defined, set it to the default value 'biased'
if (exist('biasCorrection', 'var') ~= 1)      
    biasCorrection = 'biased';   
end

% If biasCorrection doesn't start with 'b', set it to 'unbiased'
if (biasCorrection(1:1) ~= 'b')        
    biasCorrection = 'unbiased'; 
end 

% Check if the variable 'fftLength' is not defined, set it to the default value 128
if (exist('fftLength', 'var') ~= 1)      
    fftLength = 128;        
end

% Ensure that fftLength is a positive value
if (fftLength <= 0)               
    fftLength = 128;        
end 

% Check if the variable 'windowType' is not defined, set it to 0
if (exist('windowType', 'var') ~= 1)      
    windowType = 0;          
end

% Limit the maximum lag value to be within the range of 0 to sampleSize-1
maxLag = min(maxLag, sampleSize-1); 

% If fftLength is smaller than 2*maxLag+1, set it to the smallest power of 2 that is greater than or equal to sampleSize
if (fftLength < 2*maxLag + 1)   
    fftLength = 2^nextpow2(sampleSize); 
end  

%% Window Creation

% Check if the windowType is 0
if (windowType == 0)
    % If windowType is 0, create a window using Parzen window function with a length of maxLag
    window = [0; parzenwin(maxLag)];
else
     % If windowType is not 0, create a rectangular window with a length of maxLag+1
     window = ones(maxLag+1, 1);     
end

% Append zeros to the window to make its length equal to 2*maxLag+1
window = [window; zeros(maxLag, 1)];

%% Cumulants in non-redundant region

% Compute the actual overlap value based on the sample size and the given percentOverlap
actualOverlap = fix(sampleSize * percentOverlap / 100);

% Compute the advance between successive sample windows
windowAdvance  = sampleSize - actualOverlap;

% Compute the number of non-redundant records based on the total number of rows, columns, and the advance
numNonRedundantRecords  = fix((numRows * numCols - actualOverlap) / windowAdvance );

% Initialize the cumulants matrix
cumulantsMatrix = zeros(maxLag+1, maxLag+1);

% Create an index vector for selecting the samples in each window
windowIndices = (1:sampleSize)';

% Loop over each record
for recordIndex = 1:numNonRedundantRecords 
    % Extract the samples for the current window
    windowSamples = X(windowIndices);
    
    % Subtract the mean from the samples
    windowSamples = windowSamples - mean(windowSamples);
    
    % Update the index for the next window
    windowIndices = windowIndices + windowAdvance;
    
    % Loop over the lag values
    for lag = 0:maxLag
        % Compute the product of samples with lag
        laggedProducts = windowSamples(1:sampleSize-lag) .* windowSamples(lag+1:sampleSize);
        
        % Loop over the cumulant order values
        for order = lag:maxLag
            % Compute the sum of products for the current lag and order
            productSum = laggedProducts(1:sampleSize-order)' * windowSamples(order+1:sampleSize);
            
            % Apply bias correction if needed
            if (biasCorrection(1:1) == 'b')
                productSum = productSum / sampleSize;
            else
                productSum = productSum / (sampleSize - order);
            end
            
            % Update the cumulants matrix
            cumulantsMatrix(order+1, lag+1) = cumulantsMatrix(order+1, lag+1) + productSum;
        end
    end
end

% Normalize the cumulants matrix by the number of records
cumulantsMatrix = cumulantsMatrix / numNonRedundantRecords ;

%% Cumulants elsewhere by symmetry

% Apply symmetry to the cumulants matrix
symmetricCumulantsMatrix = cumulantsMatrix + tril(cumulantsMatrix, -1)';

% Extract the submatrix excluding the first row and column
cumulantsElsewhere = symmetricCumulantsMatrix(2:maxLag+1, 2:maxLag+1); 

% Initialize matrices for different cumulant components
cumulantComponent32 = zeros(maxLag, maxLag);  
cumulantComponent33 = cumulantComponent32;  
cumulantComponent34 = cumulantComponent32; 

% Loop over each order
for order = 1:maxLag
    % Extract the window samples for the current order
    windowSamples = cumulantsElsewhere(order:maxLag, order); 
    
    % Assign the window samples to the appropriate locations in the corresponding cumulant component matrices
    cumulantComponent32(maxLag+1-order, 1:maxLag+1-order) = windowSamples'; 
    cumulantComponent34(1:maxLag+1-order, maxLag+1-order) = windowSamples; 
    
    % If the order is less than maxLag, update cumulantComponent33 by considering diagonals
    if (order < maxLag) 
        windowSamples = flipud(windowSamples(2:length(windowSamples))); 
        cumulantComponent33 = cumulantComponent33 + diag(windowSamples, order) + diag(windowSamples, -order); 
    end 
end 

% Add the main diagonal elements from the first row of symmetricCumulantsMatrix to cumulantComponent33
cumulantComponent33  = cumulantComponent33 + diag(cumulantsMatrix(1, maxLag+1:-1:2)); 

% Construct the final cumulants matrix by combining cumulantComponent33, cumulantComponent32, and cumulantComponent34
cumulantsMatrix = [ [cumulantComponent33, cumulantComponent32, zeros(maxLag, 1)]; [ [cumulantComponent34; zeros(1, maxLag)] , cumulantsMatrix ] ]; 
      
%% Apply lag-domain window

% Initialize the windowed cumulants matrix
windowedCumulantsMatrix = cumulantsMatrix; 

% Check if the window type is not equal to 1 (indicating no windowing)
if (windowType ~= 1)
    % Create an index vector for the lag values
     lagIndices = (-maxLag:maxLag)';

     % Loop over each lag value
     for recordIndex=-maxLag:maxLag
         % Apply the lag-domain window to each column of the cumulants matrix
         windowedCumulantsMatrix(:,recordIndex+maxLag+1) = cumulantsMatrix(:,recordIndex+maxLag+1).* window(abs(lagIndices-recordIndex)+1) .* window(abs(lagIndices)+1)* window(abs(recordIndex)+1); 
     end 
end

%% Compute 2D-FFT, shift, and rotate for proper orientation

% Compute the 2D-FFT of the windowed cumulants matrix
bispectrum = fft2(windowedCumulantsMatrix, fftLength, fftLength); 

% Shift the bispectrum to center the axes d and r
bispectrum = fftshift(bispectrum);

% Compute the frequency axis values
if (rem(fftLength,2) == 0)
    % For even fftLength, use the appropriate frequency axis values
    waxis = (-fftLength/2:(fftLength/2-1))/fftLength;
else
    % For odd fftLength, use the appropriate frequency axis values
    waxis = (-(fftLength-1)/2:(fftLength-1)/2)/fftLength;
end 

% Plot the magnitude of the bispectrum using contour plot
contour(waxis,waxis,abs(bispectrum),4), grid on

% Set the title and labels for the plot
%title('Bispectrum estimated via the indirect method - Rectangular window')
xlabel('f1'), ylabel('f2')
end