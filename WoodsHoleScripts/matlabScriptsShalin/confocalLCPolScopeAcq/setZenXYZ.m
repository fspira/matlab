function [ status ] = setZenXYZ( tcpipZen,x,y,z )
msgdelim=';;';
timeout=5;
fwrite(tcpipZen,['-move_stage_xyz' ' ' num2str(x) ' ' num2str(y) ' ' num2str(z) msgdelim]);
status=waitforZen(tcpipZen,timeout);
%pause(delay); %500ms delay for stage to move and response to arrive.

end

