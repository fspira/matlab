function [ status ] = acquireZenImg(tcpipZen,datadir,filename)
% Acquire and save an image in datadir.
msgdelim=';;';
timeout=30*60; % Time out of 30 minutes. Tilescans take a lot of time to finish.
tobj=tic;
fwrite(tcpipZen,['-acquire_experiment "' datadir filename '"' msgdelim]);
tSendCommandToZen=toc(tobj)

tobj=tic;
status=waitforZen(tcpipZen,timeout);
tReturnStatusFromZen=toc(tobj)
end

