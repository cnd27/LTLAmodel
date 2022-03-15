function mix = getMix(p,popdat,home,nonworkmov,t,cent,movp)
%GETMIX

    %popdat.W is a LTLA by LTLA matrix of movements for work between LTLAs
    
    % Get relative sizes of populations
    Wprop = zeros(length(popdat.W),length(popdat.W));
    for i = 1:length(popdat.W)
        for j = 1:length(popdat.W)
            Wprop(i,j) = popdat.pop(i,2)/(popdat.pop(i,2)+popdat.pop(j,2));
        end
    end
    % Get proportion of each LTLA that works in each LTLA
    Ww = Wprop.*(popdat.W+popdat.W');
    Ww = Ww./repmat(sum(Ww,2),[1,length(Ww)]);
    
    for i = 1:length(Ww)
        Ww(i,i) = 1-(sum(Ww(i,:))-Ww(i,i));
    end
    % Scale away some long distance movement
    Ww2 = Ww.*exp(-cent*10^(-movp));
    
    % Scale by proportion working from home due to COVID
    Wwh = (home').*eye(length(popdat.W))+ (1-home').*Ww2;
    
    % Split by age group
    mix = zeros(3,length(popdat.W),length(popdat.W));
    if nonworkmov == 2
        for i = 1:3
            rW = repmat(p.Na(:,i),[1,length(Ww)]);
            Extra = p.Na(:,i) - sum(floor(squeeze(Wwh).*rW),2);
            Nin = floor(squeeze(Wwh).*rW);
            tmp1 = squeeze(Wwh).*rW;
            for qq = 1:length(p.Na(:,i))
                tmp = tmp1(qq,:);
                [~,ii] = maxk(tmp-floor(tmp),Extra(qq));
                Nin(qq,ii) = Nin(qq,ii) + 1;
            end
            mix(i,:,:) = Nin;
        end
    else
        if mod(t,7) == 4 || mod(t,7) == 5 % if weekend
            if nonworkmov == 0 % If no non-worker movement, no movement
                mix(1,:,:) = diag(p.Na(:,1));
                mix(2,:,:) = diag(p.Na(:,2));
                mix(3,:,:) = diag(p.Na(:,3));
            else % else random movement like Ww
                n = [popdat.pop(:,1); popdat.pop(:,2); popdat.pop(:,3)];
                rrr = rand(sum(n),1);
                n = [0; cumsum(n)];
                br = [zeros(936,1),cumsum(repmat(Ww,[3,1]),2)];
                mm = zeros(936,312);
                for i = 1:(length(n)-1)
                    mm(i,:) = histcounts(rrr((n(i)+1):n(i+1)),br(i,:));
                end
                mix(1,:,:) = mm(1:length(Ww),:);
                mix(2,:,:) = mm((length(Ww)+1):(2*length(Ww)),:);
                mix(3,:,:) = mm((2*length(Ww)+1):(3*length(Ww)),:);
            end
        else    % if not weekend
            rW = repmat(p.Na(:,2),[1,length(Ww)]);
            Extra = p.Na(:,2) - sum(floor(squeeze(Wwh).*rW),2);
            Nin = floor(squeeze(Wwh).*rW);
            tmp1 = squeeze(Wwh).*rW;
            for qq = 1:length(p.Na(:,2))
                tmp = tmp1(qq,:);
                [~,ii] = maxk(tmp-floor(tmp),Extra(qq));
                Nin(qq,ii) = Nin(qq,ii) + 1;
            end
            if nonworkmov == 0 % Only workers move (fixed)
                mix(1,:,:) = diag(p.Na(:,1));
                mix(2,:,:) = Nin;
                mix(3,:,:) = diag(p.Na(:,3));
            else % Workers move (fixed) and others random
                n = [popdat.pop(:,1); popdat.pop(:,3)];
                rrr = rand(sum(n),1);
                n = [0; cumsum(n)];
                br = [zeros(624,1),cumsum(repmat(Ww,[2,1]),2)];
                mm = zeros(624,312);
                for i = 1:(length(n)-1)
                    mm(i,:) = histcounts(rrr((n(i)+1):n(i+1)),br(i,:));
                end
                mix(1,:,:) = mm(1:length(Ww),:);
                mix(2,:,:) = Nin;
                mix(3,:,:) = mm((length(Ww)+1):(2*length(Ww)),:);
            end
        end
    end
end

