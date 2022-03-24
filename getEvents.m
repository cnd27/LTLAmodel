function events = getEvents(rates,taustep,state,nage,nvar,nexp,nvacstate,NN,neg)
%GETEVENTS Number of each event calculated

    % Get events for given rates with Poisson dist or negative binomial for infections
    if neg == 1
        pp = 0.5;
        events = cat(4,nbinrnd(rates(:,:,:,1,:)*taustep*(1-pp)/pp,1-pp),poissrnd(rates(:,:,:,2:end,:)*taustep));
        events(isnan(events)) = 0;
    else
        events = poissrnd(rates*taustep);
    end

    % Check not to many entering exposed by different variants
    if sum(sum(events(:,:,:,1,:),3) > sum(state(:,:,:,1,:),3),'all') > 0
        for i = 1:nage
            for j = 1:nvacstate
                for k = 1:NN
                    if sum(events(i,j,:,1,k),3) > sum(state(i,j,:,1,k),3)
                        vp = squeeze(events(i,j,:,1,k)./sum(events(i,j,:,1,k),3));
                        vp(isnan(vp)) = 0;
                        vp(1) = 1-sum(vp(2:end));
                        events(i,j,:,1,k) = mnrnd(sum(state(i,j,:,1,k),3),vp);
                    end
                end
            end
        end
    end
    if sum(sum(events(:,:,:,1,:),3) > sum(state(:,:,:,1,:),3),'all') > 0
        1;
    end

    events(:,:,:,2:nexp,:) = min(events(:,:,:,2:nexp,:),state(:,:,:,2:nexp,:));
    events(:,:,:,nexp+3,:) = min(events(:,:,:,nexp+3,:),state(:,:,:,nexp+2,:));
    events(:,:,:,nexp+4,:) = min(events(:,:,:,nexp+4,:),state(:,:,:,nexp+3,:));
    
    % Check not too many become infected symp/asymp
    ov = find(events(:,:,:,nexp+1,:) + events(:,:,:,nexp+2,:) > state(:,:,:,nexp+1,:));
    if ~isempty(ov)
        [i1, i2, i3, i4] = ind2sub([nage nvacstate nvar, NN], ov);% Use this
        n = zeros(length(ov),1);
        p = zeros(length(ov),1);
        for i = 1:length(ov)
            n(i) = state(i1(i),i2(i),i3(i),nexp+1,i4(i));
            p(i) = events(i1(i),i2(i),i3(i),nexp+1,i4(i))./(events(i1(i),i2(i),i3(i),nexp+1,i4(i))+events(i1(i),i2(i),i3(i),nexp+2,i4(i)));
        end
        outrnd = binornd(n,p);
        for i = 1:length(ov)
            events(i1(i),i2(i),i3(i),nexp+1,i4(i)) = outrnd(i);
            events(i1(i),i2(i),i3(i),nexp+2,i4(i)) = state(i1(i),i2(i),i3(i),nexp+1,i4(i)) - outrnd(i);
        end
    end
    if sum(isnan(events),'all') > 0
        error('Something went wrong with events')
    end

end

