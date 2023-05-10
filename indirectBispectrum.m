function Bspec = indirectBispectrum (X, nlag, nsamp,overlap, flag, nfft, wind)
%% Parameter Validation
[ly, nrecs] = size (X); 
if (ly == 1) 
    X=X(:);   
    ly = nrecs; 
    nrecs = 1;      
end 
if (exist('overlap','var') ~= 1)   
    overlap = 0;          
end 
overlap = min(99, max(overlap,0)); 
if (nrecs > 1)               
    overlap = 0;          
end 
if (exist('nsamp','var') ~= 1)     
    nsamp   = ly;         
end
if (nsamp > ly || nsamp <= 0) 
    nsamp   = ly;         
end
if (exist('flag','var') ~= 1)      
    flag    = 'biased';   
end 
if (flag(1:1) ~= 'b')        
    flag    = 'unbiased'; 
end 
if (exist('nfft','var') ~= 1)      
    nfft    = 128;        
end
if (nfft <= 0)               
    nfft    = 128;        
end 
if (exist('wind','var') ~= 1)      
    wind    = 0;          
end
nlag = min(nlag, nsamp-1); 
if (nfft  < 2*nlag+1)   
    nfft = 2^nextpow2(nsamp); 
end  

%% Window Creation
if (wind == 0) 
    window = [0; parzenwin(nlag)];
else
     window = ones(nlag+1,1);     
end 
window = [window; zeros(nlag,1)];

%% Cumulants in non-redundant region
overlap  = fix(nsamp * overlap / 100); 
nadvance = nsamp - overlap; 
nrecord  = fix ( (ly*nrecs - overlap) / nadvance ); 

c3 = zeros(nlag+1,nlag+1);
ind = (1:nsamp)';
for k=1:nrecord
    x = X(ind); x = x - mean(x);
    ind = ind + nadvance; 
    for j=0:nlag
        z = x(1:nsamp-j) .* x(j+1:nsamp); 
        for i=j:nlag
            sum = z(1:nsamp-i)' * x(i+1:nsamp); 
            if (flag(1:1) == 'b'), sum = sum/nsamp; 
            else, sum = sum / (nsamp-i); 
            end 
            c3(i+1,j+1) = c3(i+1,j+1) + sum; 
        end
    end
end
c3 = c3 / nrecord; 

%% Cumulants elsewhere by symmetry
c3 = c3 + tril(c3,-1)';
c31 = c3(2:nlag+1,2:nlag+1); 
c32 = zeros(nlag,nlag);  c33 = c32;  c34 = c32; 
for i=1:nlag 
    x = c31(i:nlag,i); 
    c32(nlag+1-i,1:nlag+1-i) = x'; 
    c34(1:nlag+1-i,nlag+1-i) = x; 
    if (i < nlag) 
        x = flipud(x(2:length(x))); 
        c33 = c33 + diag(x,i) + diag(x,-i); 
    end 
end 

c33  = c33 + diag(c3(1,nlag+1:-1:2)); 
cmat = [ [c33, c32, zeros(nlag,1)]; [ [c34; zeros(1,nlag)] , c3 ] ]; 
      
%% Apply lag-domain window
wcmat = cmat; 
if (wind ~= 1) 
     indx = (-nlag:nlag)';
     for k=-nlag:nlag
         wcmat(:,k+nlag+1) = cmat(:,k+nlag+1).* window(abs(indx-k)+1) .* window(abs(indx)+1)* window(abs(k)+1); 
     end 
end

%% Compute 2d-fft, shift and rotate for proper orientation
Bspec = fft2(wcmat, nfft, nfft); 
Bspec = fftshift(Bspec);               % axes d and r; orig at ctr

if (rem(nfft,2) == 0) 
    waxis = (-nfft/2:(nfft/2-1))/nfft;
else
    waxis = (-(nfft-1)/2:(nfft-1)/2)/nfft;
end 

contour(waxis,waxis,abs(Bspec),4), grid on 
title('Bispectrum estimated via the indirect method - Rectangular window')
xlabel('f1'), ylabel('f2')
end