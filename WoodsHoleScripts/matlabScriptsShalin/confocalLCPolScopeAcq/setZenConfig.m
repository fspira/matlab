function [ status ] = setZenConfig( tcpipZen, configName )
% Disable Zstack and time in Zen.

msgdelim=';;';
timeout=5;

fwrite(tcpipZen,['-load_config' ' ' '"' configName '"' ' ' msgdelim]);
status=waitforZen(tcpipZen,timeout);

end

