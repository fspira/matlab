function acquireZenPolStackOverTime( tcpipZen,variLC, datadir,acquisitionTime,varargin )
arg.nFrames=1;
arg.DelayBetweenT=0; % Gap between successive acquisitions.
arg.Config='';
arg.fileFormat='.lsm';
arg=parsepropval(arg,varargin{:});
nT=arg.nFrames;

setZenForPolStack(tcpipZen); %Disable z-stack, tiles, time, etc. if left on by chance.

tObj=tic; % Start a timer.
for idt=1:nT
    fprintf(1,'%s %-10s %s %-10s %s\n','Frame ', int2str(idt), 'of ', int2str(nT),':::::');
    
    I0name=['I4-0'   '_T'  num2str(idt,'%04u') '_' arg.Config arg.fileFormat];
    I135name=['I5-135'  '_T'  num2str(idt,'%04u') '_' arg.Config arg.fileFormat];
    I90name=['I6-90'  '_T'  num2str(idt,'%04u')  '_' arg.Config arg.fileFormat];
    I45name=['I7-45' '_T'  num2str(idt,'%04u') '_' arg.Config arg.fileFormat];

    %%%%%% I0 
    setLC(variLC,'I0');
    acquireZenImgFixedDelay(tcpipZen,datadir,I0name,acquisitionTime);
    t=toc(tObj);
    fprintf(1,'%s %-6s%s,',' I0', num2str(t), 's');

    %%%%%% I135
    setLC(variLC,'I135');
    acquireZenImgFixedDelay(tcpipZen,datadir,I135name,acquisitionTime);
    t=toc(tObj);            
    fprintf(1,'%s %-6s%s,',' I135', num2str(t), 's');

    %%%%%% I90
    setLC(variLC,'I90');
    acquireZenImgFixedDelay(tcpipZen,datadir,I90name,acquisitionTime);
    t=toc(tObj);            
    fprintf(1,'%s %-6s%s,',' I90', num2str(t), 's');

    %%%%%% I45
    setLC(variLC,'I45');
    acquireZenImgFixedDelay(tcpipZen,datadir,I45name,acquisitionTime);
    t=toc(tObj);            
    fprintf(1,'%s %-6s%s,',' I45', num2str(t), 's');

    pause(arg.DelayBetweenT);
end
    
end


