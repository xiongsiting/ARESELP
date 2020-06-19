function processframe(infilename,outfilename,paramfile)
% fid = fopen(paramfile);
% param = textscan(fid);
% fclose(fid);
% 
% param = param{1};


%% Step 0: Data Input
[geoinfo,echogram] = readdata(infilename);
% convert surface and bed time to row numbers 
% imDat, raw data as natural logarithm
% imAmp, which is converted from imDat to dB units
[imDat,imAmp, ysrf,ybtm] = preprocessing(geoinfo,echogram);

%% Step 1: Produce the feature image
% define the wavelet scales
scales = 3:15;
% choose the wavelet 'mexh' or 'morl'
wavelet = 'mexh';
% decide how many pixels below bed layer 
% is counted as background noise
bgSkip = 50;
peakim = peakimcwt(imAmp,scales,wavelet,ysrf,ybtm,bgSkip);

%% Step 2: Pick Seed points
% peakim : the CS or S/L ratio image 
% which is used for tracing layers.
% seed point selection
seedpt = selectseedpt(peakim);

%% Step3: layer-tracing and post-processing
% DIST = 3; BLOCKSIZE = 5;
% define parameters: distance allowance, block size, slope angle difference
DIST = 7; BLOCKSIZE = 51; SMOOTHANGLE = 90;
if nargin < 3
  params = {DIST,BLOCKSIZE,SMOOTHANGLE}; 
end

% imLayer: the image containing traced layers
% debuginfo: no. lengths, and remaining point output for debugging
tic; 
[imLayer,dinfo] = tracelayers(peakim,seedpt,params); 
toc;
tic; 
[newimLayer, labelLayer] = postprocesslayers(imLayer,DIST);
toc;

filename = split(outfilename,'.');
outtxt = [filename{1},'.txt'];
dinfo(max(labelLayer),:) = [];
dlmwrite(outtxt,dinfo);
%%

%% Step 4: Geocode layers and save them to files.
savegeolayers(outfilename,geoinfo,imAmp,newimLayer,labelLayer);

%% Steo 5: Plot layers
% load(outfilename);
% plotgeolayers(geolayers,'rand');

end
