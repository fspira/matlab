function   writePolStack(polstack,datafile,what,varargin)
% This function writes-out the polstack in uint16 format.
% The order of dimensions is assumed to be XYCT.


% Format of FluorPol stack expected by Rudolf's plugins:
% - Ratio/Anisotropy (scaled to be 12-bits).
% - Orientation (angle in degrees*10)
% - I0
% - I135
% - I90
% - I45
% - I0 (optional)
% - Total image

% Author and Copyright: Shalin Mehta, HFSP Postdoctoral Fellow
%                           Marine Biological Laboratory, Woods Hole, MA
%                           http://www.mshalin.com
% 
% License: Restricted academic use. 
% This software is free to use, but only with explicit permission by the
% author and only for academic purpose. This software constitutes a part of
% unpublished method for measuring orientation of single and ensemble of
% fluorophores.

% Ceiling is mapped to the highest value in dynamic range only during
% export.
optargs.anisoCeiling=1;
optargs.colorCeiling=1;
optargs.suffix=what;
optargs.bitDepth=16;
optargs.invertIntensity=false; % Used only for colormaps.
optargs.legend=false;
optargs.order='XYCTZ';
optargs.colorMap='hsv';
optargs.scaleIntensity=1;
optargs=parsepropval(optargs,varargin{:});

if(size(polstack,3)~=7 && size(polstack,3)~=8 )
    error('3rd dimension of Polstack must equal 7 or 8.');
end
X=size(polstack,2);
Y=size(polstack,1);
C=size(polstack,3);
Z=size(polstack,4);
T=size(polstack,5);

[pathstr, name, ext]=fileparts(datafile);
name(name == ' ')='_';
polStackName=[pathstr '/' name   optargs.suffix  ext];
orientfileName=[pathstr '/' name  optargs.suffix ext];

switch(what)
    case ('OME-TIFF')
        % OME-TIFF is much (>10x) slower than the saveastiff call below.
        
        polstack(:,:,1,:)=((2^16-1)/optargs.anisoCeiling)*polstack(:,:,1,:);
        polstack(:,:,2,:)=100*(180/pi)*polstack(:,:,2,:);
        bfsave(uint16(polstack),polStackName,optargs.order);
    case('PolStack')
        switch(optargs.bitDepth)
            case 8
                polstack(:,:,1,:)=((2^8-1)/optargs.anisoCeiling)*polstack(:,:,1,:);
                polstack(:,:,2,:)=(180/pi)*polstack(:,:,2,:);
                %polstack=reshape(polstack,[Y X C*Z*T]);
%                 metadata_text = char([
%                 ['ImageJ=1.48d', 10],...
%                 ['images=', num2str(C*Z*T), 10],...
%                 ['channels=', num2str(C), 10],...                
%                 ['slices=', num2str(Z), 10],...
%                 ['frames=', num2str(T), 10],...
%                 ['hyperstack=true', 10],...
%                 ['mode=grayscale', 10],...                
%                 ['loop=false', 10]]);     
%                 options.ImageDescription=metadata_text;
%                % options.comp='lzw';
%                 saveastiff(uint8(polstack),polStackName,options);
                polstack=reshape(polstack,size(polstack,1),size(polstack,2),size(polstack,3),size(polstack,4),1); % Reshape according to XYCZT convention.
                exportHyperStack(uint8(polstack),polStackName);
                %writeHyperStackTag(polStackName,C,Z,T);              
            case 10
                
            case 12
                polstack(:,:,1,:)=((2^12-1)/optargs.anisoCeiling)*polstack(:,:,1,:);
                polstack(:,:,2,:)=10*(180/pi)*polstack(:,:,2,:);
                %polstack=reshape(polstack,[Y X C*Z*T]);
                %saveastiff(uint16(polstack),polStackName);
                polstack=reshape(polstack,size(polstack,1),size(polstack,2),size(polstack,3),size(polstack,4),1); % Reshape according to XYCZT convention.
                exportHyperStack(uint16(polstack),polStackName);
            case 14
                
            case 16
                polstack(:,:,1,:)=((2^16-1)/optargs.anisoCeiling)*polstack(:,:,1,:);
                polstack(:,:,2,:)=100*(180/pi)*polstack(:,:,2,:);
                %polstack=reshape(polstack,[Y X C*Z*T]);
                %saveastiff(uint16(polstack),polStackName);
                polstack=reshape(polstack,size(polstack,1),size(polstack,2),size(polstack,3),size(polstack,4),1); % Reshape according to XYCZT convention.
                exportHyperStack(uint16(polstack),polStackName);

                %writeHyperStackTag(polStackName,C,Z,T);              
               
            otherwise
                error('Bit-depth can be 8,10,12,14 or 16.');
        end
        % Following code is >10x slower on Windows machine than the above
        % code.
%         for idT=1:size(polstack,4)
%             aniso=uint16((2^16-1)*polstack(:,:,1,idT));
% 
%             if(idT == 1) % If this is the very first file of the polstack create a new file.
%                 imwrite(aniso,datapath);
%             else
%                 imwrite(aniso,datapath,'WriteMode','append');
%             end
% 
%             azim=uint16(100*(180/pi)*polstack(:,:,2,idT));
%             imwrite(azim,datapath,'WriteMode','append');
% 
%             avg=uint16(polstack(:,:,3,idT));
%             imwrite(avg,datapath,'WriteMode','append');
% 
%             I0=uint16(polstack(:,:,4,idT));
%             imwrite(I0,datapath,'WriteMode','append');
% 
%             I135=uint16(polstack(:,:,5,idT));
%             imwrite(I135,datapath,'WriteMode','append');
% 
%             I90=uint16(polstack(:,:,6,idT));
%             imwrite(I90,datapath,'WriteMode','append');
% 
%             I45=uint16(polstack(:,:,7,idT));
%             imwrite(I45,datapath,'WriteMode','append');
%        end
        case('OrientationMap')
        % Map azimuth to hue (periodic at boundaries 0 and 1).
        % Map saturation to 1 (could map anisotropy to saturation).
        % Map value (brightness to total intensity).
            aniso=squeeze(polstack(:,:,1,:));
            % Clip the anisotropy to colorCeiling.
            aniso(aniso>optargs.colorCeiling)=optargs.colorCeiling;
            azim=squeeze(polstack(:,:,2,:));
            avg=optargs.scaleIntensity*squeeze(polstack(:,:,3,:));
            avg=avg/(2^optargs.bitDepth-1); % Normalize w.r.t. the dynamic range.
            if(optargs.invertIntensity) %Useful for color visualization of diattenuation data.
                avg=max(avg(:))-avg;
            end
            OrientationMap=pol2color(aniso,azim,avg,optargs.colorMap,'legend',optargs.legend);
            options.color=true;
            saveastiff(uint8(255*OrientationMap),orientfileName,options);
            
end


end

function writeHyperStackTag(filename,channelN,slicesN,framesN)
% Writing following Metadata entry to TIFF stack converts it into
% Hyperstack. The default write order in hyperstack is XYCZT.
% This function contributed by Amitabh Verma.

imagesN=channelN*slicesN*framesN;
% http://metarabbit.wordpress.com/2014/04/30/building-imagej-hyperstacks-from-python/
metadata_text = char([
                ['ImageJ=1.48o', 10],...
                ['images=', num2str(imagesN), 10],...
                ['channels=', num2str(channelN), 10],...                
                ['slices=', num2str(slicesN), 10],...
                ['frames=', num2str(framesN), 10],...
                ['hyperstack=true', 10],...
                ['mode=grayscale', 10],...                
                ['loop=false', 10]]);
        
    t = Tiff(fullfile(filename), 'r+');
    % Modify the value of a tag.
    % http://www.digitalpreservation.gov/formats/content/tiff_tags.shtml
    % t.setTag('Software',['OI-DIC',' ',num2str(version_number),'x']);
    % t.setTag('Model', 'Lumenera Infinity 3M');
    
    % set tag
    t.setTag('ImageDescription', metadata_text);
    
    % write tag to file
    t.rewriteDirectory();
    
    % close tiff
    t.close()
end