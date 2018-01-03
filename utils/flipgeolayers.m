function geolayers = flipgeolayers(geolayers)

    lonStart = geolayers.longitude(1); 
    lonEnd = geolayers.longitude(end);
    latStart = geolayers.latitude(1);
    latEnd = geolayers.latitude(end);
    
    if lonStart > lonEnd || (lonStart == lonEnd && latStart > latEnd)
        fields = fieldnames(geolayers);
        for i = 1:size(fields,1)
            if isvector(geolayers.(fields{i})) && ...
                    size(geolayers.(fields{i}),2) == geolayers.num_trace...
                    && ~isstruct(geolayers.(fields{i}))
                geolayers.(fields{i}) = fliplr(geolayers.(fields{i}));
            end
        end
        
        for j = 1:geolayers.num_layer
            fields = fieldnames(geolayers.layer(j));
            for i = 1:size(fields,1)
                if isvector(geolayers.layer(j).(fields{i})) && ...
                        size(geolayers.layer(j).(fields{i}),2) > 1 && ...
                    ~isstruct(geolayers.layer(j).(fields{i}))
                    geolayers.layer(j).(fields{i}) = ...
                        fliplr(geolayers.layer(j).(fields{i}));
                end
            end
        end
    end

end