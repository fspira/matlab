function processConfocalPolTimeSeries( directory,Frames,varargin )
arg.Channel=1; % Select the channel to process.
%arg.Position=1;
arg.registerPol=false; % Do pol-channels need registration?
arg.bin=1;
arg.normFactors=[1 1 1 1];
arg.suffix=''; %Suffix after the slice name.
arg.prefix='';
arg.BlackLevel=0;
arg.outputFormat=NaN; %Bit-depth has no effect on how data is read. It only affects the format in which computed results are written.
arg.displayStatus=true;
arg=parsepropval(arg,varargin{:});


% Process one time-point at a time and export the computed results. Avoid
% accumulation in RAM.
filenames={};

for idF=Frames
%    textprogressbar(['Frame#' num2str(idF,'%u') ' in' directory ':']);
% Bioformats reader keeps printing filenames and indicates progress.
    
        [I0,I45,I90,I135]=getPolChannels(directory,idF,arg.Channel,arg.prefix); % Change this function as naming convention changes.
         
        if isnan(arg.outputFormat)
            outputFormat=class(I0);
        else
            outputFormat=arg.outputFormat;
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
        ComputeFluorAnisotropy(I0,I45,I90,I135,...
        'anisotropy','BlackLevel',arg.BlackLevel,'normFactors',arg.normFactors);
    
        [anisoPath,orientPath,avgPath]=writeComputedChannels(aniso,orient,avg,directory,idF,arg.Channel,arg.suffix,outputFormat);
        filenames=cat(1,filenames,{anisoPath,orientPath,avgPath}');
        % The order of the files ends up being channels, slices, and
        % frames.This stack can be re-ordered after import.
        
%        textprogressbar(round(idZ/max(arg.Zslices)));
%    textprogressbar(['DONE Frame#' num2str(idF,'%u')]);
    
end

% Write out the list so that Fiji can import the results as a stack.
listi = fopen([directory '/PolOutput.txt'],'w');
for idN=1:numel(filenames)
    fprintf(listi,'%s\n',filenames{idN}); 
end
fclose(listi);

end

function [I0,I45,I90,I135]=getPolChannels(directory,Frame,Channel,prefix)  % Change this function as naming convention changes.
    
    I0T=(Frame-1)*4+1;
    I45T=(Frame-1)*4+2;
    I90T=(Frame-1)*4+3;
    I135T=(Frame-1)*4+4;
    
    I0file=[directory prefix 't' num2str(I0T,'%02u') '.lsm'];
    I45file=[directory prefix 't' num2str(I45T,'%02u') '.lsm'];
    I90file=[directory  prefix 't' num2str(I90T,'%02u') '.lsm'];
    I135file=[directory  prefix 't' num2str(I135T,'%02u') '.lsm'];

    r0=bfGetReader(I0file);
    I0=bfGetPlane(r0,Channel); 

    r45=bfGetReader(I45file);
    I45=bfGetPlane(r45,Channel); 

    r90=bfGetReader(I90file);
    I90=bfGetPlane(r90,Channel); 

    r135=bfGetReader(I135file);
    I135=bfGetPlane(r135,Channel); 
    
end

function  [anisoPath,orientPath,avgPath]=writeComputedChannels(aniso,orient,avg,directory,Frame,Channel,suffix,bitDepth)
Tstr=num2str(Frame,'%03u'); Chstr=num2str(Channel,'%u');
anisoPath=[directory '/I1-aniso'  '_T'  Tstr '_Ch' Chstr suffix '.tif'];
orientPath=[directory '/I2-orient'  '_T'  Tstr '_Ch' Chstr suffix '.tif'];
avgPath=[directory '/I3-avg'  '_T'  Tstr '_Ch' Chstr  suffix '.tif'];

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