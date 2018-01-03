% function peakim = peakimcwt(im,scales,wavelet,ysrf,ybtm,app,bgSkip)
function peakim = peakimcwt(im,scales,wavelet,ysrf,ybtm,bgSkip)
%--------------------------------------------------------------------------
% Generate the feat image by using average of CWT coefficient
%
% Parameters
% ----------
% im : 2-D image
%   MCoRDS echogram
% ysrf : 1 * n array
%   Row numbers of surface echoes
% ybtm : 1 * n array
%   Row numbers of bottom echoes
% wavelet : string
%   Name of mother wavelet, 'mexh' or 'mexh'
% scales : 1 * n array
%   Withs of the mother wavelet, [1,4,8,16]
% bgSkip : int
%   Pixels below bed reflectors to start calculating background noise, 50
% 
% See also
% --------
% read_data.m
% preprocessing.m
%
% Returns
% -------
% peakim : 2-D image
%   Feature image as same dimension as 'im'
%   Pixel values of peaks are wavelet coefficient
%   Non-peaks are indicated by zero.
%
% Examples
% --------
% [geoinfo,echogram] = read_data(inputfilename);
% [im_nl,im_db, ysrf,ybtm] = preprocessing(geoinfo,echogram);
% peakim = find_peaks_cwt(im,ysrf,ybtm);
%--------------------------------------------------------------------------

    %% Define parameters
    if nargin < 6
        bgSkip = 50; % set to use app method or not
    end
    
%     if nargin < 6
%         app = 0;
%     end
%     
    if nargin < 5
        ybtm = ones(1,size(im,2))*size(im,1);
    end
    
    if nargin < 4
        ysrf = ones(1,size(im,2));      
    end

    if nargin < 3
        wavelet = 'mexh';
    end
    
    if nargin < 2
        scales = 1:16;
    end

%     coefim = zeros(size(im));
    peakim = zeros(size(im));
    
%     if app == 1
%         kurtoim = zeros(size(im));
%         ratioim = zeros(size(im));
%     end
    %% Start processing
    for i = 1: size(im,2)
%         if app == 1
%             [peakim(:,i),~,kurtoim(:,i),ratioim(:,i)] = findpeakscwt(im(:,i),...
%                 scales,wavelet,...
%                 ysrf(i),ybtm(i),bgSkip);
%             
%         else
            peakim(:,i) = findpeakscwt(im(:,i), scales,wavelet,...
                ysrf(i),ybtm(i),bgSkip);
%         end
    end
    
%     if app == 1
%                     
%         % linear fit of ratio (x) and Kurtosis (y)
%         % y = p1 * x + p2;
%         [y,x] = size(ratioim);
%         allratio = reshape(ratioim,[x*y,1]);
%         allkurto = reshape(kurtoim,[x*y,1]);
%         allratio(allratio == 0) = [];
%         allkurto(allkurto == 0) = [];
%         pft = polyfit(allkurto,allratio,1);
%         peakim(ratioim < pft(1) * kurtoim + pft(2)) = 0;
%     end
    
%     peakim = scale2minmax(peakim);
%     peakim_nan = peakim;
%     peakim_nan(peakim == 0) = nan;
%     minval = min(peakim_nan(:));
%     maxval = max(peakim_nan(:));
%     peakim = (peakim_nan - minval)/(maxval-minval);
%     peakim(isnan(peakim)) = 0;

end


