function acquireZenPolStack( tcpipZen,variLC, x,y,z, datadir,DelayBetweenXY, DelayBetweenZ,varargin )
nZ=size(z,2);
nP=size(z,1);
tObj=tic; % STart a timer.
arg.PolConfig=true;
arg.Configs='ch1';
arg.fileFormat='.lsm';
arg.delayBetweenChannels=5; % Delay in seconds when changing channels.
arg.nFrames=1;
arg.DelayBetweenT=0; % Gap between successive acquisitions.
arg=parsepropval(arg,varargin{:});

nC=numel(arg.Configs);
nT=arg.nFrames;
% Acquisition order is assumed to be Z, Channels , Position, Time.
% Switching channels takes a lot of time on confocal and can lead to errors
% in acquisition.

for idt=1:nT
        fprintf(1,'%s %-10s %s %-10s %s\n','Frame ', int2str(idt), 'of ', int2str(nT),':::::');

    for idp=1:nP
        fprintf(1,'%s %-10s %s %-10s %s\n','Tile ', int2str(idp), 'of ', int2str(nP),':::::');
        posdir=[datadir '/P' num2str(idp,'%03u') '/'];
     
        for idc=1:nC % For each channel
                setZenConfig(tcpipZen,arg.Configs{idc});
                fprintf(1,'%s %-5s %s \n','Channel ', int2str(idc), arg.Configs{idc});
                pause(arg.delayBetweenChannels); % Let Zen settle in new configuration.
                setZenForPolStack(tcpipZen); %Disable z-stack, tiles, time, etc. if left on by chance.
            for idz=1:nZ
                % Move stage.
                zstatus=setZenXYZ( tcpipZen,x(idp),y(idp),z(idp,idz));        
                if(zstatus)
                    fprintf(1,'%s %-10s %s %-10s %s','Z ', int2str(idz), 'of ', int2str(nZ), ':');
                    if(arg.PolConfig(idc)) % Check if this config is polarization sensitive or not.
                        I0name=['I4-0' '_Z' num2str(idz,'%.3d')  '_T'  num2str(idt,'%04u') '_' arg.Configs{idc} arg.fileFormat];
                        I135name=['I5-135' '_Z' num2str(idz,'%.3d')  '_T'  num2str(idt,'%04u') '_' arg.Configs{idc} arg.fileFormat];
                        I90name=['I6-90' '_Z' num2str(idz,'%.3d')  '_T'  num2str(idt,'%04u')  '_' arg.Configs{idc} arg.fileFormat];
                        I45name=['I7-45' '_Z' num2str(idz,'%.3d')  '_T'  num2str(idt,'%04u') '_' arg.Configs{idc} arg.fileFormat];
                        I0bleachname=['I8-0' '_Z' num2str(idz,'%.3d')  '_T'  num2str(idt,'%04u') '_' arg.Configs{idc} arg.fileFormat];
                        %%%%%% I0 
                        setLC(variLC,'I0');
                        acquireZenImg(tcpipZen,posdir,I0name);
                        t=toc(tObj);
                        fprintf(1,'%s %-6s%s,',' I0', num2str(t), 's');

                        %%%%%% I135
                        setLC(variLC,'I135');
                        acquireZenImg(tcpipZen,posdir,I135name);
                        t=toc(tObj);            
                        fprintf(1,'%s %-6s%s,',' I135', num2str(t), 's');

                        %%%%%% I90
                        setLC(variLC,'I90');
                        acquireZenImg(tcpipZen,posdir,I90name);
                        t=toc(tObj);            
                        fprintf(1,'%s %-6s%s,',' I90', num2str(t), 's');

                        %%%%%% I45
                        setLC(variLC,'I45');
                        acquireZenImg(tcpipZen,posdir,I45name);
                        t=toc(tObj);            
                        fprintf(1,'%s %-6s%s,',' I45', num2str(t), 's');

                        %%%%%% I0bleach/ Circular illumination.
                        setLC(variLC,'I0bleach');
                        acquireZenImg(tcpipZen,posdir,I0bleachname);
                        t=toc(tObj);            
                        fprintf(1,'%s %-6s%s\n',' I0bleach', num2str(t),'s');

                    else
                        Iname=['I' '_Z' num2str(idz,'%.3d') '_' arg.Configs{idc} arg.fileFormat];
                        acquireZenImg(tcpipZen,posdir,Iname);
                    end
    
                else
                    error(['Couldn''t move the stage to position: x=' num2str(x) ',y=' num2str(y) ',z=' num2str(z(idp,idz))]);
                end

                pause(DelayBetweenZ);

            end
        end
            pause(DelayBetweenXY);

    end
    
pause(arg. DelayBetweenT);
end

end

