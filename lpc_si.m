function [si lpc_map] = lpc_si(im, C, Beta_k, scales, norient, B)

if (~exist('C'))
   C = 2;
end

if (~exist('Beta_k'))
   Beta_k = 1e-4;
end

if (~exist('scales'))
scales   = [1 3/2 2]; w = [1 -3 2];
% scales = [1 2 3]; w = [1 -4 3];
% scales = [1 2 4]; w = [1 -3 2];
end

if (~exist('norient'))
    norient  = 8;
end

nscale   = length(scales);
[row col]= size(im);

if (~exist('B'))
    B = round(min(row, col)/16);
end

im = double(im);
filter  = logGabor_2D(im,norient,nscale,scales,0.33);
imfft   = fft2(im);
s_lpc   = ones(row, col, norient);

for o = 1:norient
    for s = 1:nscale
        M(:,:,s) = ifft2(imfft .* filter{s,o}); 
        s_lpc(:,:,o) = s_lpc(:,:,o).*(M(:,:,s).^w(s));
    end
    e = abs(M(:,:,1));
    e_center = e(B+1:end-B, B+1:end-B);
    e_mean = mean(e_center(:));
    e_std = std(e_center(:));
    T = e_mean + 2*e_std;
    e = max(0, e-T);
    energy(:,:,o) = e;
end
s_lpc_map = cos(angle(s_lpc));
s_lpc_map(s_lpc_map<0) = 0;
lpc_map = (sum(s_lpc_map.*energy,3))./(sum(energy,3) + C);
lpc_map_center = lpc_map(B+1:end-B, B+1:end-B);

sorted_si = sort(lpc_map_center(:),'descend')';
N = length(sorted_si);
u = exp(-((0:(N-1))/(N-1))/Beta_k); 
si = sum(sorted_si.*u)/sum(u);

return;
