function Bspec = directBispectrum (X, nfft, wind, nsamp, overlap)
%% Parameter Validation
[ly, nrecs] = size(X);
if (ly == 1) 
    X = X(:);  
    ly = nrecs; 
    nrecs = 1; 
end
if (exist('nfft','var') ~= 1)            
    nfft = 128; 
end
if (exist('overlap','var') ~= 1)      
    overlap = 50;  
end
overlap = min(99,max(overlap,0));
if (nrecs > 1)                  
    overlap =  0;  
end
if (exist('nsamp','var') ~= 1)          
    nsamp = 0;   
end
if (nrecs > 1)                    
    nsamp = ly;  
end
if (nrecs == 1 && nsamp <= 0)
nsamp = fix(ly/ (8 - 7 * overlap/100));
end
if (nfft  < nsamp)   
    nfft = 2^nextpow2(nsamp); 
end

overlap  = fix(nsamp * overlap / 100);
nadvance = nsamp - overlap;
nrecs    = fix ( (ly*nrecs - overlap) / nadvance);


%% 2D Window Creation
if (exist('wind', 'var') ~= 1) 
    wind = 5; 
end
[m,n] = size(wind);
window = wind;

if (max(m,n) == 1)     % scalar: wind is size of Rao-Gabr window
    winsize = wind;
    if (winsize < 0) 
        winsize = 5; 
    end        % the window length L
    winsize = winsize - rem(winsize,2) + 1;  % make it odd
    if (winsize > 1)
        mwind   = fix (nfft/winsize);            % the scale parameter M
        lby2    = (winsize - 1)/2;

        theta  = -lby2:lby2;
        opwind = ones(winsize,1) * (theta .^2);       % w(m,n)=m^2
        opwind = opwind + opwind' + theta' * theta;   % m^2 + n^2 + mn
        opwind = 1 - (2*mwind/nfft)^2 * opwind;       %
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
Bspec    = zeros(nfft,nfft);

mask = hankel(1:nfft,[nfft,1:nfft-1] );   % the hankel mask (faster)
locseg = (1:nsamp)';
for krec = 1:nrecs
    xseg   = X(locseg);
    Xf     = fft(xseg-mean(xseg), nfft)/nsamp;
    CXf    = conj(Xf);
    Bspec  = Bspec + (Xf * Xf.') .*reshape(CXf(mask), nfft, nfft);
    locseg = locseg + nadvance;
end

Bspec = fftshift(Bspec)/(nrecs);

%% Frequency-Domain Smoothing
if (winsize > 1)
    lby2 = (winsize-1)/2;
    Bspec = conv2(Bspec,opwind);
    Bspec = Bspec(lby2+1:lby2+nfft,lby2+1:lby2+nfft);
end


%% Contour Plot of Magnitude Bispectum
if (rem(nfft,2) == 0)
    waxis = (-nfft/2:(nfft/2-1))'/nfft;
else
    waxis = (-(nfft-1)/2:(nfft-1)/2)'/nfft;
end

contour(waxis,waxis,abs(Bspec),4),grid on 
title('Bispectrum estimated via the direct (FFT) method')
xlabel('f1'), ylabel('f2')
end