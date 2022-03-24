function [state] = doVaccine(state,p,t,vac,Nexp,Nvar,Nage,Nvacstates)
%DOVACCINE Moves people to correct vaccine class

    % Get number to be vaccinated at time t for each dose (Age, LTLA)
    v1 = vac(:,:,t,1); 
    v2 = vac(:,:,t,2);
    
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
