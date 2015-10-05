function processConfocalPolData( directory,varargin )
% processConfocalPolData( directory,varargin ) iterates over dimensions and
% outputs computed anisotropy, average and orientation. Since TIFF file
% sizes are limited to ~4GB, output for each slice is written in separate
% file.
arg.Zslices=1; 
arg.Frames=1;
arg.Channel=1; % Select the channel to process.
%arg.Position=1;
arg.registerPol=false; % Do pol-channels need registration?
arg.bin=1;
arg.normFactors=[1 1 1 1];
arg.suffix=''; %Suffix after the slice name.
arg.prefix=''; %Prefix before the time-stamp for data acquired by polling when files are written by Zen.
arg.BlackLevel=90;
arg.bitDepth=NaN; %Bit-depth has no effect on how data is read. It only affects the format in which computed results are written.
arg.displayStatus=true;
arg.acqMethod='zen-controller';
arg=parsepropval(arg,varargin{:});


% Process one slice at a time and export the computed results. Avoid
% accumulation in RAM.
filenames={};

for idF=arg.Frames
%    textprogressbar(['Frame#' num2str(idF,'%u') ' in' directory ':']);
% Bioformats reader keeps printing filenames and indicates progress.
    

    for idZ=arg.Zslices
        
        % Change this decision tree and corresponding funciton for
        % different methods of acqusition and naming conventions.
        switch(arg.acqMethod)
            case {'polling','Polling','POLLING'}
             [I0,I45,I90,I135]=getPolChannelsPolling(directory,idF,arg.Channel,arg.prefix,max(arg.Frames));
            otherwise
            [I0,I45,I90,I135]=getPolChannels(directory,idZ,idF,arg.Channel,arg.suffix);
        end
         
        if isnan(arg.bitDepth)
            bitDepth=class(I0);
        else
            bitDepth=arg.bitDepth;
        end
        
        if(arg.bin~=1)
            scale=1/arg.bin;
            I0=imresize(I0,scale,'bilinear');
            I45=imresize(I45,scale,'bilinear');
            I90=imresize(I90,scale,'bilinear');
            I135=imresize(I135,scale,'bilinear');
       % I0bleach=imresize(I0bleach,scale,'bilinear');
        end
        
        if(arg.registerPol)
            X=size(I0,2); Y=size(I0,1); 
            I45(:,:,idz)=imregphasecor(I0(:,:,idz),I45(:,:,idz),1:X,1:Y,'translation');
            I90(:,:,idz)=imregphasecor(I0(:,:,idz),I90(:,:,idz),1:X,1:Y,'translation');
            I135(:,:,idz)=imregphasecor(I0(:,:,idz),I135(:,:,idz),1:X,1:Y,'translation');
            %I0bleach(:,:,idz)=imregphasecor(I0(:,:,idz),I0bleach(:,:,idz),1:X,1:Y,'translation');        
        end
       
        [orient, aniso, avg]=...
        ComputeFluorAnisotropy_1(I0,I45,I90,I135,...
        'anisotropy','BlackLevel',arg.BlackLevel,'normFactors',arg.normFactors);
    
        [anisoPath,orientPath,avgPath]=writeComputedChannels(aniso,orient,avg,directory,idZ,idF,arg.Channel,arg.suffix,bitDepth);
        filenames=cat(1,filenames,{anisoPath,orientPath,avgPath}');
        % The order of the files ends up being channels, slices, and
        % frames.This stack can be re-ordered after import.
        
%        textprogressbar(round(idZ/max(arg.Zslices)));
    end
%    textprogressbar(['DONE Frame#' num2str(idF,'%u')]);
    
end

% Write out the list so that Fiji can import the results as a stack.
listi = fopen([directory '/analysis/PolOutput.txt'],'w');
for idN=1:numel(filenames)
    fprintf(listi,'%s\n',filenames{idN}); 
end
fclose(listi);

end

function [I0,I45,I90,I135]=getPolChannels(directory,Slice,Frame,Channel,suffix)
    
    I0file=[directory '/I4-0' '_Z' num2str(Slice,'%03u') '_T'  num2str(Frame,'%04u') suffix '.lsm'];
    I45file=[directory '/I7-45' '_Z' num2str(Slice,'%03u') '_T'  num2str(Frame,'%04u') suffix '.lsm'];
    I90file=[directory '/I6-90' '_Z' num2str(Slice,'%03u') '_T'  num2str(Frame,'%04u') suffix '.lsm'];
    I135file=[directory '/I5-135' '_Z' num2str(Slice,'%03u') '_T'  num2str(Frame,'%04u') suffix '.lsm'];

    r0=bfGetReader(I0file);
    I0=bfGetPlane(r0,Channel); 

    r45=bfGetReader(I45file);
    I45=bfGetPlane(r45,Channel); 

    r90=bfGetReader(I90file);
    I90=bfGetPlane(r90,Channel); 

    r135=bfGetReader(I135file);
    I135=bfGetPlane(r135,Channel); 
    
end


% Copy and modify this function to read pol-channels per frame.
function [I0,I45,I90,I135]=getPolChannelsPolling(directory,Frame,Channel,prefix,totalFrames)  % Change this function as naming convention changes.

% Construct file-names
I0T=(Frame-1)*4+1;
I45T=(Frame-1)*4+2;
I90T=(Frame-1)*4+3;
I135T=(Frame-1)*4+4;


getDirectory =  dir(char(directory));
lengthDirectory = length(getDirectory);

clear listDirectory
clear sortedList

for subrun = 1:length(getDirectory)-4
    laufName = getDirectory(subrun+4).name;
    listDirectory{subrun} =  laufName;
    suffixStrSort = findstr( laufName,'.lsm');
    prefixStrSort = findstr( laufName,'T');
    
    
    sortArray(subrun) = str2num(laufName(prefixStrSort+1:suffixStrSort-1));
end

[sortOut sortIdx] = sort(sortArray);

sortedList{1} = getDirectory(4).name;

for subrun = 1:length(sortIdx)
    sortedList{subrun+1} = listDirectory{sortIdx(subrun)};
end

I0file = [char(directory) '\' char(sortedList(I0T))];
I45file=[char(directory) '\' char(sortedList(I45T))];
I90file=[char(directory) '\' char(sortedList(I90T))];
I135file=[char(directory) '\' char(sortedList(I135T))];



%  if(totalFrames<10)
%      I0file=[directory prefix 't' num2str(I0T,'%1u') '.lsm'];
%      I45file=[directory prefix 't' num2str(I45T,'%1u') '.lsm'];
%      I90file=[directory  prefix 't' num2str(I90T,'%1u') '.lsm'];
%      I135file=[directory  prefix 't' num2str(I135T,'%1u') '.lsm'];
%  elseif(totalFrames<100)
%      I0file=[directory prefix 't' num2str(I0T,'%02u') '.lsm'];
%      I45file=[directory prefix 't' num2str(I45T,'%02u') '.lsm'];
%      I90file=[directory  prefix 't' num2str(I90T,'%02u') '.lsm'];
%      I135file=[directory  prefix 't' num2str(I135T,'%02u') '.lsm'];
%  elseif(totalFrames<1000)
%      I0file=[directory prefix 't' num2str(I0T,'%03u') '.lsm'];
%      I45file=[directory prefix 't' num2str(I45T,'%03u') '.lsm'];
%      I90file=[directory  prefix 't' num2str(I90T,'%03u') '.lsm'];
%      I135file=[directory  prefix 't' num2str(I135T,'%03u') '.lsm'];
%  else
%      I0file=[directory prefix 't' num2str(I0T,'%04u') '.lsm'];
%      I45file=[directory prefix 't' num2str(I45T,'%04u') '.lsm'];
%      I90file=[directory  prefix 't' num2str(I90T,'%04u') '.lsm'];
%      I135file=[directory  prefix 't' num2str(I135T,'%04u') '.lsm'];
%  end


% Use bio-formats to read pixels.
r0=bfGetReader(I0file);
I0=bfGetPlane(r0,Channel);

r45=bfGetReader(I45file);
I45=bfGetPlane(r45,Channel);

r90=bfGetReader(I90file);
I90=bfGetPlane(r90,Channel);

r135=bfGetReader(I135file);
I135=bfGetPlane(r135,Channel);

end

function  [anisoPath,orientPath,avgPath]=writeComputedChannels(aniso,orient,avg,directory,Slice,Frame,Channel,suffix,bitDepth)
Zstr=num2str(Slice,'%03u'); Tstr=num2str(Frame,'%04u'); Chstr=num2str(Channel,'%u');
anisoPath=[directory '/analysis/I1-aniso' '_Z' Zstr '_T'  Tstr '_Ch' Chstr suffix '.tif'];
orientPath=[directory '/analysis/I2-orient' '_Z' Zstr '_T'  Tstr '_Ch' Chstr suffix '.tif'];
avgPath=[directory '/analysis/I3-avg' '_Z' Zstr '_T'  Tstr '_Ch' Chstr  suffix '.tif'];

    switch(bitDepth)
        case 'uint8'
            aniso=(2^8-1)*aniso;
            orient=(180/pi)*orient;
            imwrite(uint8(aniso),anisoPath);
            imwrite(uint8(orient),orientPath);
            imwrite(uint8(avg),avgPath);
            
        case 'uint16'
            aniso=(2^16-1)*aniso;
            orient=100*(180/pi)*orient;
            imwrite(uint16(aniso),anisoPath);
            imwrite(uint16(orient),orientPath);
            imwrite(uint16(avg),avgPath);
            
        case 'single'
            saveastiff(anisoPath,single(aniso));
            saveastiff(orientPath*(180/pi),single(orient));
            saveastiff(avgPath,single(avg));
    end
end