function ICs = getICs(p,Nage,Nvacstates,Nvar,Nexp,vac,vacdays,infe,rec,tstart,alpw,alpn,impdistalpha,newvar)

%GETICS Summary of this function goes here
%   Detailed explanation goes here

    ICs = zeros(Nage,Nvacstates,Nvar,4+Nexp,size(p.Na,1));
    % Get initial vaccianted
    if tstart > 0
        vacd = sum(vac(:,:,1:(tstart-1),:),3);
    else
        vacd = zeros(Nage,size(p.Na,1),1,2);
    end
    vacd = min(vacd,repmat(p.Na',[1,1,1,2]));
    
    % Set numbers in vaccinated classes
    ICs(:,3,1,1,:) = permute(vacd(:,:,:,2),[1,3,4,5,2]);
    ICs(:,2,1,1,:) = permute(vacd(:,:,:,1)-vacd(:,:,:,2),[1,3,4,5,2]);
    ICs(:,1,1,1,:) = permute(p.Na,[2,3,4,5,1]) - sum(ICs(:,:,1,1,:),2);
    
    % Set numbers in recovered classes
    recsym = round(permute((p.Na').*repmat(rec,[1,312]),[1,3,4,5,2]).*permute(vacdays,[1,3,4,5,2]));
    recsym = min(recsym,ICs(:,:,1,1,:));
    ICs(:,:,1,(Nexp+4),:) = recsym(:,:,:,1,:);
    ICs(:,:,1,1,:) = ICs(:,:,1,1,:) - ICs(:,:,1,Nexp+4,:);

    % Load infection distribution October 2020
    load('./data/LineListDataAge.mat')
    tmp = squeeze(casesLTLAwa(39,:,:));
    tmp(tmp==0) = 1;
    locinf = (sum(p.N)*infe*tmp/sum(sum(tmp)));
    infprob = repmat([repmat((1/p.epsilon)/(Nexp/p.epsilon+1/p.gamma),[3,Nexp]), (p.d(:,1)/p.gamma)/(Nexp/p.epsilon+1/p.gamma), ((1-p.d(:,1))/p.gamma)/(Nexp/p.epsilon+1/p.gamma)],[1,1,p.Nltlas]);
    infs = round(permute(locinf,[2,3,1]).*infprob);

    % Set numbers in infected classes
    ICs(:,1,1,2:(Nexp+3),:) = infs;
    ICs(:,1,1,1,:) = ICs(:,1,1,1,:) - (sum(ICs(:,1,1,(2:Nexp+3),:),4));
    
    % Seed new variant
    if tstart == 279
        if impdistalpha == 1
            if newvar == 1
                load('./data/vardists.mat')
                impdistalp = alpdist;
            else   
                load('./data/impdistalp.mat')
            end
            wwt = histcounts(randsample(1:312,alpn,1,impdistalp),0.5:1:312.5);
            place = find(wwt>0);
            number = wwt(wwt>0);
            symp = binornd(reshape(number,[length(number),1]),p.d(2,1));
            ICs(2,1,2,Nexp+2,place) = symp;
            ICs(2,1,2,Nexp+3,place) = reshape(number,[length(number),1]) - symp;
            ICs(2,1,1,1,place) = ICs(2,1,1,1,place) - sum(ICs(2,1,2,(Nexp+2):(Nexp+3),place),4);
        else
            l2 = load('./data/regionnum.mat');
            regionnum = l2.regionnum;
            ww = zeros(312,1);
            ww(regionnum==1|regionnum==3) = alpw;
            ww(~(regionnum==1|regionnum==3)) = 1;
            wwt = histcounts(randsample(1:312,alpn,1,p.N/sum(p.N).*ww),0.5:1:312.5);
            place = find(wwt>0);
            number = wwt(wwt>0);
            symp = binornd(reshape(number,[length(number),1]),p.d(2,1));
            ICs(2,1,2,Nexp+2,place) = symp;
            ICs(2,1,2,Nexp+3,place) = reshape(number,[length(number),1]) - symp;
            ICs(2,1,1,1,place) = ICs(2,1,1,1,place) - sum(ICs(2,1,2,(Nexp+2):(Nexp+3),place),4);
        end
    end
end

