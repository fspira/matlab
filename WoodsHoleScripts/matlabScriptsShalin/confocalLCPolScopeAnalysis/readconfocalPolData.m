function [ I0,I45,I90,I135, I0bleach, black ] = readconfocalPolData( directory,varargin )

optargs.Zslices=1;
optargs.Time=1;
optargs.Channel=1; % For multi-channel data, we need to select one channel.
optargs.Position=1;
optargs.registerPol=false; % Do pol-channels need registration?
optargs.bin=1;
optargs=parsepropval(optargs,varargin{:});


I0=getZStack(directory, 'I4-0',optargs.Zslices,optargs.Channel,optargs.Time,optargs.Position);
I45=getZStack(directory, 'I7-45',optargs.Zslices,optargs.Channel,optargs.Time,optargs.Position);
I90=getZStack(directory, 'I6-90',optargs.Zslices,optargs.Channel,optargs.Time,optargs.Position);
I135=getZStack(directory, 'I5-135',optargs.Zslices,optargs.Channel,optargs.Time,optargs.Position);
I0bleach=getZStack(directory, 'I8-0',optargs.Zslices,optargs.Channel,optargs.Time,optargs.Position);

if(optargs.bin~=1)
    scale=1/optargs.bin;
    I0=imresize(I0,scale,'bilinear');
    I45=imresize(I45,scale,'bilinear');
    I90=imresize(I90,scale,'bilinear');
    I135=imresize(I135,scale,'bilinear');
    I0bleach=imresize(I0bleach,scale,'bilinear');
end

 X=size(I0,2); Y=size(I0,1); Z=size(I0,3);


if(optargs.registerPol)
    for idz=1:Z
        I45(:,:,idz)=imregphasecor(I0(:,:,idz),I45(:,:,idz),1:X,1:Y,'translation');
        I90(:,:,idz)=imregphasecor(I0(:,:,idz),I90(:,:,idz),1:X,1:Y,'translation');
        I135(:,:,idz)=imregphasecor(I0(:,:,idz),I135(:,:,idz),1:X,1:Y,'translation');
        I0bleach(:,:,idz)=imregphasecor(I0(:,:,idz),I0bleach(:,:,idz),1:X,1:Y,'translation');        
    end
end


I0=reshape(I0,[Y X 1 Z]);
I45=reshape(I45,[Y X 1 Z]);
I90=reshape(I90,[Y X 1 Z]);
I135=reshape(I135,[Y X 1 Z]);
I0bleach=reshape(I0bleach,[Y X 1 Z]);


black=min([min(I0(:)) min(I45(:)) min(I90(:)) min(I135(:)) min(I0bleach(:))]);

end

function Zstack=getZStack(directory,PolChannel,Zslices,Channel,Time,Position)
    
    % Try new format. October 1, 2014 onwards.
    searchstr=[directory '/' PolChannel '_Z*' '_T'  num2str(1,'%04u') '_'  '*.lsm'];
    files=dir(searchstr); % New format. 
    if(isempty(files)) % Try old format.
        searchstr=[directory '/' PolChannel '_Z*' '_P'  num2str(Position,'%.3d') '*.lsm'];
        files=dir(searchstr);
    end
    
    if(isempty(files))
        error(['No file matches this path:' searchstr]);
    end
    
    filesRead=files(Zslices); % sub-sample the file names if necessary.

    
    r=bfGetReader([directory filesRead(1).name]); % Using the reader avoids the overhead, when the file contains multiple channels/tiles.
%    r.setSeries(Tile-1);  % Needed only if LSM file contains multiple
%    tiles.
    Zstack=double(bfGetPlane(r,Channel*Time));
    Zstack=repmat(Zstack,[1 1 numel(Zslices)]);
   
    % Read other slices.
    for idz=2:numel(Zslices)
        r=bfGetReader([directory filesRead(idz).name]);
%    r.setSeries(Tile-1);  % Needed only if LSM file contains multiple
%    tiles.
        Zstack(:,:,idz)=double(bfGetPlane(r,Channel*Time));     
    end
end

% function voxels=getVoxelsFromFile(file,Channel)
% [~,~,extension]=fileparts(file);
% switch(extension(2:end)) %Discard the dot.
%     case {'tif','TIF','TIFF','tiff'}
%         TIFFobj=TIFFStack(file);
%         voxels=TIFFobj(:,:,Channel); %Assume third dimension is channel.
%     case {'lsm','LSM','czi','CZI'}
%         datawithmeta=bfOpen3DVolume(file);
%         voxels=datawithmeta{1}{1}(:,:,Channel); %Assume third dimension is channel.
% end
% end
