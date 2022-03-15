function vacdays = getVacDays(tstart,vac,p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    vacdaysn = zeros(3,312,size(vac,3),3);

    cumvac = cumsum(vac,3);
    cumvac = min(cumvac,permute(p.Na,[2,1]));

    vacdaysn(:,:,:,3) = cumvac(:,:,:,2);
    vacdaysn(:,:,:,2) = cumvac(:,:,:,1)-cumvac(:,:,:,2);
    vacdaysn(:,:,:,1) = permute(p.Na,[2,1]) - (cumvac(:,:,:,1));

    vacdays = squeeze(sum(vacdaysn(:,:,1:(tstart-1),:),3))./sum(squeeze(sum(vacdaysn(:,:,1:(tstart-1),:),3)),3);

end

