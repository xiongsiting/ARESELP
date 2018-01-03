function [imDat,imAmp, ysrf,ybtm] = preprocessing(geoinfo,echogram,trimrows)
%--------------------------------------------------------------
% Calcuate echogram in dB units and the row numbers of surface and bottom layers
%
% Parameters
% ----------
% geoinfo : struct 
%   Containing geoinformation
% echogram : 2-D image
%   MCoRDS echogram
%  
% Returns
% -------
% imDat : 2-D image
%   raw data as natural logarithm
% imAmp : 2-D image
%   which is converted from imDat to dB units
% ysrf : 1 * n array
%   row numbers of surface layer
% ybtm : 1 * n array
%   row numbers of bottom layer
% 
% See also
% --------
%   read_data.m
%
% Examples
% --------
% [geoinfo, echogram] = read_data('Data_20110329_02_019.mat');
% [imDat,imAmp, ysrf,ybtm] = preprocessing(geoinfo,echogram)
%---------------------------------------------------------------

    if nargin < 3
        trimrows = 300; %trimrows = 300, 0 for pre-icebridge
    end
    % convert the format of the echogram
    if sum(echogram(:)) < 0
        % pre *** the MCoRDS data has already converted to dB
        imDat = [];
        imAmp = echogram;
    else
        imDat = echogram;
        imAmp = 10*log10(echogram);
    end
    % trim the bottom weird area of the echogram
    imAmp(size(imAmp,1) - trimrows + 1:end,:) = [];
    if ~isempty(imDat), imDat(size(imDat,1) - trimrows + 1:end,:) = [];end
%-------------------------------------------------    
    %% Extract the surface and bottom layers
    Surface = geoinfo.traveltime_surface;
    Bottom = geoinfo.traveltime_bottom;
    Time = geoinfo.time_range;
    ysrf = zeros(size(Surface)); % row number of Surface along track
    ybtm = zeros(size(Bottom));  % row number of Bottom along track

    for i = 1:length(Surface)
        idxSurf = find(Time >= Surface(i));
        ysrf(i) = idxSurf(1);
        idxBtm = find(Time <= Bottom(i));
        if isempty(idxBtm)
            ybtm(i) = 1080;
        else
            ybtm(i) = idxBtm(end);
        end
    end
    % ysrf and ybtm contains the row numbers of surface and bottom echoes    
end