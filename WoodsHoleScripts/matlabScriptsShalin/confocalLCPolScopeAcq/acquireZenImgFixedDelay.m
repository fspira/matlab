function [ status ] = acquireZenImgFixedDelay(tcpipZen,datadir,filename,waittime)
% Acquire and save an image in datadir.
msgdelim=';;';
%tobj=tic;
fwrite(tcpipZen,['-acquire_experiment "' datadir filename '"' msgdelim]);
%tSendCommandToZen=toc(tobj)

%tobj=tic;
pause(waittime);
%tReturnStatusFromZen=toc(tobj)
end

