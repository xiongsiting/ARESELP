function batchconnectframes(datalist)

[nums,inlist,outlist] = loadfilelist(datalist);
% for one segment
for s = 8:12
    tic;
    segments = inlist(s,:);
    basenames = cell(nums(s),1);
    filelist = cell(nums(s),1);

    % processing each frame
    for j = 1:nums(s)
        
        fname = segments{j};
        tmp = split(fname,'.');
        filelist{j} = ['LAYERS-',fname];

        if j == 1, fistsgbasename = tmp{1};end
            
        if j == nums(s)
            lastframeno = split(tmp{1},'_');
            lastframeno = lastframeno{4};
            outfilename = ['LAYERS-',fistsgbasename,'-',lastframeno,'.mat'];
        end
        
        processframe(fname,filelist{j});
    end

    if ~isempty(outfilename)
        save(outfilename,'filelist');
    else
        return;
    end

    frame_layers_stack = cell(size(filelist,1),1);

    for i = 1:length(filelist)
        frame_layers_stack{i} = load(filelist{i});
        if i == 1
            frame1 = frame_layers_stack{i}; 
            geolayers = frame1.geolayers;
        end
        addfrm = frame_layers_stack{i}; 
        geolayers = connectframes(geolayers,addfrm.geolayers);
    end
    save(outfilename,'geolayers','frame_layers_stack','-append');
%     plotgeolayers(geolayers,'rand');
    toc;
end

end
