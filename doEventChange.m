function state = doEventChange(state,events,Nexp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    state(:,:,1,1,:) = state(:,:,1,1,:) - sum(events(:,:,:,1,:),3);
    state(:,:,:,2,:) = state(:,:,:,2,:) + events(:,:,:,1,:);
    state(:,:,:,Nexp+1,:) = state(:,:,:,Nexp+1,:) - events(:,:,:,Nexp+1,:);
    state(:,:,:,Nexp+2,:) = state(:,:,:,Nexp+2,:) + events(:,:,:,Nexp+1,:);
    state(:,:,:,Nexp+1,:) = state(:,:,:,Nexp+1,:) - events(:,:,:,Nexp+2,:);
    state(:,:,:,Nexp+3,:) = state(:,:,:,Nexp+3,:) + events(:,:,:,Nexp+2,:);
    state(:,:,:,Nexp+2,:) = state(:,:,:,Nexp+2,:) - events(:,:,:,Nexp+3,:);
    state(:,:,:,Nexp+4,:) = state(:,:,:,Nexp+4,:) + events(:,:,:,Nexp+3,:);
    state(:,:,:,Nexp+3,:) = state(:,:,:,Nexp+3,:) - events(:,:,:,Nexp+4,:);
    state(:,:,:,Nexp+4,:) = state(:,:,:,Nexp+4,:) + events(:,:,:,Nexp+4,:);

end

