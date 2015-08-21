function [ status ] = saveZenConfig( tcpipZen, configName )
% Disable Zstack and time in Zen.

msgdelim=';;';
timeout=5;

fwrite(tcpipZen,['-save_config' ' ' '"' configName '"' ' ' msgdelim]);
status=waitforZen(tcpipZen,timeout);

end

