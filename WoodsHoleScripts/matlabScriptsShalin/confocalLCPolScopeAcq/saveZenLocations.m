function [ status ] = saveZenLocations( tcpipZen, filePath )
% Disable Zstack and time in Zen.

msgdelim=';;';
timeout=5;

fwrite(tcpipZen,['-export_marked_locations' ' ' '"' filePath '"' ' ' msgdelim]);
status=waitforZen(tcpipZen,timeout);

end

