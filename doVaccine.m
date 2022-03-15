function [state,surgestate] = doVaccine(state,p,t,vac,Nexp,Nvar,Nage,Nvacstates,surge,surgeboost,surgethresh,delay,surgestate,omitdd,maxper,svar,threshtype)
%DOVACCINE 

    % Get number to be vaccinated at time t for each dose (Age, LTLA)
    v1 = vac(:,:,t,1); 
    v2 = vac(:,:,t,2);
    if surge == 1
        vv1 = zeros(size(v1,1),size(v1,2));
        vv2 = zeros(size(v2,1),size(v2,2));
        if threshtype == 1
            surgingnow = squeeze(sum(sum(sum(state(:,:,svar,(Nexp+2):(Nexp+3),:),1),2),4)./sum(sum(sum(sum(state(:,:,:,(Nexp+2):(Nexp+3),:),1),2),3),4));
        elseif threshtype == 2
            surgingnow = squeeze(sum(sum(sum(sum(state(:,:,:,(Nexp+2):(Nexp+3),:),1),2),3),4)./p.N);
        elseif threshtype == 3
            load('./data/regionnum.mat')
            surgingnow = regionnum == 6;
        end
        tosurge = surgingnow > surgethresh;
        
        % Maximum surging
        maxsurge = round(maxper*312);
        [~,tmpi] = sort(surgingnow,'descend');
        currentsurge = (surgestate(:,t)>0);
        currentsurgen = sum(surgestate(:,t)>0);
        currentsurgeh = (surgestate(:,t)==1|surgestate(:,t)==2);
        currentsurgehn = sum(surgestate(:,t)==1|surgestate(:,t)==2);
        tosurgen = zeros(312,1);
        tosurgen(currentsurgeh) = 1;
        lefttosurge = min(maxsurge+omitdd- currentsurgehn,sum(logical(tosurge)&~logical(tosurgen)));
        tmpis = setdiff(tmpi,find(currentsurgeh),'stable');
        tosurgen(tmpis((1+omitdd):(lefttosurge))) = 1;
        tosurgen(tmpis(1:min(omitdd,sum(tosurge)))) = 2;

        tosurge = tosurgen;
        
        if delay > 0
            surgestate(tosurge==1 & surgestate(:,t)==0 & (surgestate(:,t-1)==2|surgestate(:,t-1)==3),t) = 3;
            surgestate(tosurge==1 & surgestate(:,t)==0,(t+delay):(t+13+delay)) = 2;
            surgestate(tosurge==1 & surgestate(:,t)==0,t:(t+delay-1)) = 1;
            surgestate(tosurge==2,t) = -1;
        else
            surgestate(tosurge==1 & surgestate(:,t)==0 & (surgestate(:,t-1)==2|surgestate(:,t-1)==3),t) = 3;
            surgestate(tosurge==1 & surgestate(:,t)==0,t:(t+13)) = 2;
            surgestate(tosurge==2,t) = -1;
        end

        surgestateV = [zeros(312,14) surgestate];
        
        vv1(:,surgestateV(:,t)>0) = min(round(v1(:,surgestateV(:,t)>0)*surgeboost),squeeze(sum(sum(state(:,1,:,[1 5],surgestateV(:,t)>0),3),4)));
        if sum(sum(vv1)) > sum(sum(v1))
            lowersurge = fminsearch(@(x)abs(sum(sum(v1))-sum(sum(round(vv1(:,surgestateV(:,t)>0)*x)))),1);
            vv1(:,surgestateV(:,t)>0) = round(vv1(:,surgestateV(:,t)>0)*lowersurge);    
        end
        vv1(:,surgestateV(:,t)<1) = v1(:,surgestateV(:,t)<1);
        if sum(sum(vv1)) > sum(sum(v1))
            lowersurge2 = fminsearch(@(x)abs(sum(sum(v1))-sum(sum(vv1(:,surgestateV(:,t)>0)))-sum(sum(round(vv1(:,surgestateV(:,t)<1)*x)))),1);
            vv1(:,surgestateV(:,t)<1) = round(vv1(:,surgestateV(:,t)<1)*lowersurge2);
        end
        vv1(:,surgestateV(:,t)==1) = v1(:,surgestateV(:,t)==1);
        
        vv2(:,surgestateV(:,t)>0) = min(round(v2(:,surgestateV(:,t)>0)*surgeboost),squeeze(sum(sum(state(:,2,:,[1 5],surgestateV(:,t)>0),3),4)));
        if sum(sum(vv2)) > sum(sum(v2))
            lowersurge3 = fminsearch(@(x)abs(sum(sum(v2))-sum(sum(round(vv2(:,surgestateV(:,t)>0)*x)))),1);
            vv2(:,surgestateV(:,t)>0) = round(vv2(:,surgestateV(:,t)>0)*lowersurge3);
        end
        vv2(:,surgestateV(:,t)<1) = v2(:,surgestateV(:,t)<1);
        if sum(sum(vv2)) > sum(sum(v2))
            lowersurge4 = fminsearch(@(x)abs(sum(sum(v2))-sum(sum(vv2(:,surgestateV(:,t)>0)))-sum(sum(round(vv2(:,surgestateV(:,t)<1)*x)))),1);
            vv2(:,surgestateV(:,t)<1) = round(vv2(:,surgestateV(:,t)<1)*lowersurge4);
        end
        vv2(:,surgestateV(:,t)==1) = v2(:,surgestateV(:,t)==1);

        v1 = max(0,vv1);
        v2 = max(0,vv2);
                
    end
    
    if sum(v1+v2,'all') > 0
        % Check not too many vaccinated
        v1  = min(v1,squeeze(sum(state(:,1,:,1,:)+state(:,1,:,Nexp+4,:),3)));
        v2  = min(v2,squeeze(sum(state(:,2,:,1,:)+state(:,2,:,Nexp+4,:),3)));
        
        % Number of vaccinations in each (Age, LTLA, Variant, sus recI or recA, dose)
        vEventsN = zeros(Nage,2,Nvar,2,size(v1,2));
        v1 = v1(:);
        v2 = v2(:);
        
        % All random numbers
        rrr = rand(sum(sum(v1+v2)),1);
        n1 = [0; cumsum(v1)];
        n2 = [0; cumsum(v2)];
        for i = 1:length(v1)
            % Probability of S, RS, RA for each variant for dose 1 in each (age, LTLA)
            if Nvar == 3
                br1tmp = [0 state(mod(i-1,3)+1,1,1,1,ceil(i/3)) state(mod(i-1,3)+1,1,2,1,ceil(i/3)) state(mod(i-1,3)+1,1,3,1,ceil(i/3)) state(mod(i-1,3)+1,1,1,Nexp+4,ceil(i/3)) state(mod(i-1,3)+1,1,2,Nexp+4,ceil(i/3)) state(mod(i-1,3)+1,1,3,Nexp+4,ceil(i/3))];
            elseif Nvar == 2
                br1tmp = [0 state(mod(i-1,3)+1,1,1,1,ceil(i/3)) state(mod(i-1,3)+1,1,2,1,ceil(i/3)) state(mod(i-1,3)+1,1,1,Nexp+4,ceil(i/3)) state(mod(i-1,3)+1,1,2,Nexp+4,ceil(i/3))];
            end
            br1 = cumsum(br1tmp/sum(br1tmp));
            br1(isnan(br1))=0;
            a1 = histcounts(rrr((n1(i)+1):n1(i+1)),br1);
            a1=min(br1tmp(2:end),a1);
            % Probability of S, RS, RA for each variant for dose 2 in each (age, LTLA)
            if Nvar == 3
                br2tmp = [0 state(mod(i-1,3)+1,2,1,1,ceil(i/3)) state(mod(i-1,3)+1,2,2,1,ceil(i/3)) state(mod(i-1,3)+1,2,3,1,ceil(i/3)) state(mod(i-1,3)+1,2,1,Nexp+4,ceil(i/3)) state(mod(i-1,3)+1,2,2,Nexp+4,ceil(i/3)) state(mod(i-1,3)+1,2,3,Nexp+4,ceil(i/3))];
            elseif Nvar == 2
                br2tmp = [0 state(mod(i-1,3)+1,2,1,1,ceil(i/3)) state(mod(i-1,3)+1,2,2,1,ceil(i/3)) state(mod(i-1,3)+1,2,1,Nexp+4,ceil(i/3)) state(mod(i-1,3)+1,2,2,Nexp+4,ceil(i/3))];
            end
            br2 = cumsum(br2tmp/sum(br2tmp));
            br2(isnan(br2))=0;
            a2 = histcounts(rrr((n2(i)+1):n2(i+1)),br2);
            a2=min(br2tmp(2:end),a2);
            vEventsN(mod(i-1,3)+1,1,:,:,ceil(i/3)) = reshape(a1,[1,1,Nvar,2]);
            vEventsN(mod(i-1,3)+1,2,:,:,ceil(i/3)) = reshape(a2,[1,1,Nvar,2]);
        end
        % Do vaccination
        for k = 1:2
            state(:,k,:,1,:) = state(:,k,:,1,:) - vEventsN(:,k,:,1,:);
            state(:,k+1,:,1,:) = state(:,k+1,:,1,:) + vEventsN(:,k,:,1,:);
            state(:,k,:,Nexp+4,:) = state(:,k,:,Nexp+4,:) - vEventsN(:,k,:,2,:);
            state(:,k+1,:,Nexp+4,:) = state(:,k+1,:,Nexp+4,:) + vEventsN(:,k,:,2,:);
        end
           
    end
end
