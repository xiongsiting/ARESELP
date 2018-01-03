function geolayers = savegeolayers(outfilename,geoinfo,echogram,...
    echogramLayer,lbLayer,withim)

    if nargin < 6
        withim = 'withim';
    end
    C = 3e8; epsilon = 3.15;

    nlayer = length(lbLayer);
    geoinfo.num_layer = 0;
          
    for i = 1:nlayer
        [y,x] = find(echogramLayer == lbLayer(i));
        if isempty(x)
            continue;
        end
        geoinfo.num_layer = geoinfo.num_layer + 1;
        indices = sub2ind(size(echogramLayer),y,x);

        depth = ones(1,geoinfo.num_trace)*nan;
        echo_intensity = depth;
        elevation = depth;
%         if size(geoinfo.traveltime_surface,1) < size(geoinfo.traveltime_surface,2)
%             depth(x) = (geoinfo.time_range(y) - geoinfo.traveltime_surface(x)')*C/2/sqrt(epsilon);
%         else
        depth(x) = (geoinfo.time_range(y) - geoinfo.traveltime_surface(x))*C/2/sqrt(epsilon);
%         end
        try
        elevation(x) = geoinfo.elevation_surface(x) - depth(x);
        catch
            x;
        end
        echo_intensity(x) = echogram(indices);
        geoinfo.layer(geoinfo.num_layer).depth = depth;
        geoinfo.layer(geoinfo.num_layer).elevation = elevation;
        geoinfo.layer(geoinfo.num_layer).echo_intensity = echo_intensity;
    end
    geolayers = geoinfo;
    
    if strcmp(withim,'withim')
        geolayers.echogram = echogram;
    end
    % Uncomment these for flip the layers according to East-West,
    % South-North order
    % This is not good for further connect frames
    % so,instead, we realise this flipping when plotting in plotlayers.m
    % geolayers = flipgeolayers(geolayers);
    
    if ~strcmp(outfilename,'')
        save(outfilename,'geolayers');
        
    end
end

