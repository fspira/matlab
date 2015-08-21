function [ status ] = setZenLocations( tcpipZen, filepath )

msgdelim=';;';
timeout=5;

fwrite(tcpipZen,['-import_marked_locations' ' ' '"' filepath '"' ' ' msgdelim]);
status=waitforZen(tcpipZen,timeout);

end

