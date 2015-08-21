function exportConfocalPolData(directory,varargin)
% exportConfocalPolData(directory,params,value).
% params                values
% 'outputType'          'OrientationMap','PolStack','Vectors'
% 'otherChannels'       true/false. Include other channels in PolStack.
% 'Polraw'              true/false. Include raw images in PolStack.
% 'Polcomputed'         true/false. Include computed results in PolStack.

arg.outputType='OrientationMap';
arg.Channel=1;
arg.PolOutput=true;
arg.NonPol=false; % If set yes, all other channels are included.
arg.suffixNonPol='';
arg.PolRaw=false;
arg.Zslices=1;
arg.Frames=1;
arg.lengthPropToavg=true;
arg.lengthPropToaniso=true;
arg.avgCeiling=NaN; % Compute default value for avgCeiling based on bitDepth.
arg.anisoCeiling=NaN; % Since we do not have access to whole stack, use half the dynamic range for ceiling. 
arg.inputFormat=NaN;
arg.suffix='';
arg=parsepropval(arg,varargin{:});

 if(isnan(arg.avgCeiling) || isnan(arg.anisoCeiling))
    error('Since this function processes one slice at a time, it cannot estimate appropriate ceilings. They must be supplied as inputs.');
 end

 % Prepare to write the files during the loop.
[path,dirname]=fileparts(directory);
outputprefix=[path '/' dirname '_' arg.outputType];
options.append=true;

% Before starting to write, delete the existing files.
switch(arg.outputType)
    case 'OrientationMap'
        options.color=true;
        orientationFile=[outputprefix '.tif'];
        if exist(orientationFile,'file')
            delete(orientationFile);
        end
        %varargout{1}=orientationFile;
    case 'PolStack'
        options.color=false;
        polstackFile=[outputprefix '.tif'];
        if exist(polstackFile,'file')
            delete(polstackFile);
        end

        %varargout{1}=hyperstackFile;
    case 'Vectors'
        options.color=false;
        xvecFile=[outputprefix 'X.tif'];
        if exist(xvecFile,'file')
            delete(xvecFile);
        end
        
        yvecfile=[outputprefix 'Y.tif'];
        if exist(yvecfile,'file')
            delete(yvecfile);
        end
        
        anisoFile=[outputprefix 'aniso.tif'];
        if exist(anisoFile,'file')
            delete(anisoFile);
        end
        
        avgFile=[outputprefix   'intensity.tif'];
        if exist(avgFile,'file')
            delete(avgFile);
        end
    
    otherwise %Assume raw  images.
        outputfile=[outputprefix 'Ch' int2str(arg.Channel) '.tif'];
        delete(outputfile);
end


for idF=arg.Frames
%    textprogressbar(['Frame#' num2str(idF,'%u') ' in' directory ':']);
% Bioformats reader keeps printing filenames and indicates progress.
    for idz=arg.Zslices

        [aniso,orient,avg]=getPolOutputs(directory,idz,idF,arg.Channel,arg.suffix);
        if isnan(arg.inputFormat) % Determine the format during the first call.
            arg.inputFormat=class(aniso);
        end

        switch(arg.outputType)
         case 'PolStack'
                % ImageJ hyperstack assumes 'XYCZT' sequence. Therefore,
                % our appending files in CZT order is correct.
                % Order is other channels, pol-outputs, followed by raw pol
                % channels.
                channels=[];
                
                if(arg.NonPol)
                    otherChannels=getOtherChannels(directory,idz,idF,arg.suffixNonPol);
                    channels=cat(3,channels,otherChannels);
                end
                
                if(arg.PolOutput)
                    channels=cat(3,channels,aniso,orient,avg);
                end
                
                if(arg.PolRaw)
                    raw=getPolRaw(directory,idz,idF,arg.Channel,arg.suffix);
                    channels=cat(3,channels,raw);
                end
                saveastiff(channels,polstackFile,options);                    
            case 'OrientationMap'
                 % Convert to double
                  aniso=single(aniso); orient=single(orient); avg=single(avg);
                switch(arg.inputFormat)
                    case 'uint8'
                        aniso=aniso/(2^8-1);
                        orient=(pi/180)*orient;
                    case 'uint16'
                        aniso=aniso/(2^16-1);
                        orient=(pi/180)*(orient/100);
                    case 'single'
                        orient=(pi/180)*orient;
                end                
                orientationmap=pol2color(aniso,orient,avg,'sbm','anisoCeiling',arg.anisoCeiling,'avgCeiling',arg.avgCeiling,'legend',true);
                saveastiff(uint8(255*orientationmap),orientationFile,options);
           
            case 'Projections'
                
            case 'Vectors'
                switch(arg.inputFormat)
                    case 'uint8'
                        aniso=aniso/(2^8-1);
                        orient=(pi/180)*orient;
                    case 'unit16'
                        aniso=aniso/(2^16-1);
                        orient=(pi/180)*(orient/100);
                    case 'single'
                        orient=(pi/180)*orient;
                end
                [Xvector,Yvector,anisoVector,avgVector]=pol2vectors(single(aniso),single(orient),single(avg),'avgCeilig',arg.avgCeiling,'anisoCeiling',arg.anisoCeiling);                

        end
    end

end

end

function [aniso,orient,avg]=getPolOutputs(directory,Z,F,Ch,suffix)
    Tstr=num2str(F,'%04u'); Chstr=num2str(Ch,'%u'); Zstr=num2str(Z,'%03u');
    aniso=imread([directory '/' 'analysis' '/' 'I1-aniso' '_Z' Zstr '_T' Tstr '_Ch' Chstr suffix '.tif']);
    orient=imread([directory '/' 'analysis' '/' 'I2-orient' '_Z' Zstr '_T' Tstr '_Ch' Chstr suffix '.tif']);
    avg=imread([directory '/' 'analysis' '/' 'I3-avg' '_Z' Zstr '_T' Tstr '_Ch' Chstr suffix '.tif']);
end

function channels=getOtherChannels(directory,Z,F,suffix)
    Tstr=num2str(F,'%04u');  Zstr=num2str(Z,'%03u');
    channels=bfopen([directory '/' 'I' '_Z' Zstr '_T' Tstr suffix '.lsm']);
end