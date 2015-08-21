function [ status ] = setZenForPolStack( tcpipZen )
% Disable Zstack and time in Zen.

msgdelim=';;';
timeout=5;

fwrite(tcpipZen,['-set_experiment_actions' ' ' num2str(0) ' ' num2str(0) msgdelim]);
status=waitforZen(tcpipZen,timeout);

end

