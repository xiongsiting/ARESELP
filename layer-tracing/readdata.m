function [geoinfo,echogram] = readdata(filename)
%--------------------------------------------------------------
% Read geoinfo and other Meta info from MATFILE
%
% Parameters
% ----------
% filename : mat file
%   MCoRDS Mat file
%
% See also
% --------
%   netcdf_to_mat.m
%   convert2mat.m
%  
% Returns
% -------
% geoinfo : struct 
%   Containing geoinformation
% echogram : 2-D image
%   MCoRDS echogram
%
% Examples
% --------
% [geoinfo, echogram] = read_data('Data_20110329_02_019.mat');
%---------------------------------------------------------------
    
    C = 3e8; EPSILON = 3.15;
    % extract geoinformation and other Meta data from datasets
    load(filename);
    echogram = Data;
    
    fdnames = {'distance','elevation_bed','elevation_surface',...
        'latitude','longitude','num_layer','num_trace','thickness',...
        'time_gps','traveltime_surface','x','y','layer',...
        'traveltime_bottom','time_range'}; %last two fields added by XST
    
    ntcfield = 11;
    ntrace = size(Data,2);
    traceData = zeros(ntcfield,ntrace);
    
    traceData(4,:) = Latitude;
    traceData(5,:) = Longitude;
    [x,y] = polarstereo_fwd(Latitude,Longitude,...
        6378137.0,0.08181919,70,-45);
%     distance = sqrt(x.^2 + y.^2);
%     traceData(1,:) = distance;
    traceData(9,:) = x;
    traceData(10,:) = y;
    traceData(7,:) = GPS_time;
    traceData(8,:) = Surface;
    traceData(11,:) = Bottom;
    traceData(6,:) = (Bottom - Surface)*C/2/sqrt(EPSILON); % thickness
    traceData(3,:) = Elevation - Surface*C/2; % elevation_surface
    traceData(2,:) = traceData(3,:) - traceData(6,:);
    
    if traceData(1,1) == 0
        geoinfo = struct(fdnames{1},[]);
    else
        geoinfo = struct(fdnames{1},traceData(1,:));
    end
    
    ndata = 2;
    for i = 2:length(fdnames)
        switch i
            case 6
                geoinfo.num_layer = [];
%                 nfd = nfd + 1;
            case 7
                geoinfo.num_trace = ntrace;
%                 nfd = nfd + 1;
            case 13
                geoinfo.layer = [];
%                 nfd = nfd + 1;
            case 15
                geoinfo.time_range = Time;
            otherwise
                geoinfo.(fdnames{i}) = traceData(ndata,:);
                ndata = ndata + 1;
        end
        
    end

end