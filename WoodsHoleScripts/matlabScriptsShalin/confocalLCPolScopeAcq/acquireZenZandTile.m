function acquireZenZandTile( tcpipZen, x,y,z, datadir,DelayBetweenXY, DelayBetweenZ )
%  acquireZenZandTile( tcpipZen, x,y,z, datadir,DelayBetweenXY, DelayBetweenZ )
% tcpipZen: tcpip pipe to Zen-controller.
% x,y: vectors of same size specifying position of the stage.
% z: vector specifying focus range.
% datadir: directory in which each image is saved individually.
% DelayBetweenZ: Delay between moving focus. 
% DelayBetweenXY: Delay between moving stage.
% Above delays are needed if writing a large file to the disk.

nZ=numel(z);
nP=numel(x);
tObj=tic; % STart a timer.
for idp=1:nP
    fprintf(1,'%s %-10s %s %-10s %s\n','X= ', num2str(x(idp)), 'Y= ', num2str(y(idp)),':::::');

    for idz=1:nZ
        % Move stage.
        zstatus=setZenXYZ( tcpipZen,x(idp),y(idp),z(idz));        
        if(zstatus)
            fprintf(1,'%s %-10s %s','Z= ', num2str(z(idz)), ':');
            acquireZenImg(tcpipZen,datadir,['I' '_Z' num2str(idz,'%.3d') '_P' num2str(idp,'%.3d')]);
            t=toc(tObj);
            fprintf(1,'%-6s%s,', num2str(t), 's');

        else
            error(['Couldn''t move the stage to position: x=' num2str(x) ',y=' num2str(y) ',z=' num2str(z(idz))]);
        end
        
        pause(DelayBetweenZ);

    end
        pause(DelayBetweenXY);

end

end

