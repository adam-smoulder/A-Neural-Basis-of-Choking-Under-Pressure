function [ ] = parsave( filename, data )
% save for a single variable to be used in parallel computing shenanigans
t = getCurrentTask();
save([filename num2str(mod(t.ID,2)+1)],'data');
end

