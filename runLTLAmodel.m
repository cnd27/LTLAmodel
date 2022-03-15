function [state, cases, casesI,surgestate] = runLTLAmodel(p,vac,Nvar,Nexp,Nvacstates,Nage,nonworkmov,tstart,tend,prev,prec3,mob,agemix,da,tau,R0,tr,alpw,alpn,deln,mobscale,impdistalpha,deltastart,neg,newvar,ltlaprecmob,surge,surgeboost,surgethresh,delay,omitdd,maxper,svar,threshtype)

    % Set parameters
    p.R0 = R0;
    p.tau = tau;
    p.da = da; 
    p.tr = tr;
    p.d = p.d.*(da');
    ageprop = p.Na./sum(p.Na,2);

    bbb = getBetaInit(p,agemix,ageprop);
    p.beta = zeros(p.Nltlas,Nvar);
    p.beta(:,1) = bbb;
    p.beta(:,2) = tr(1)*p.beta(:,1);
    p.beta(:,3) = tr(1)*tr(2)*p.beta(:,1);
    
    % Get times 
    tspan = [tstart tend];
    tvals = sort([tspan(1):tspan(2),(tspan(1)+0.7):(tspan(2)-0.3)]);
    tvalsi = tspan(1):tspan(2);
    step = [p.delta 1-p.delta]; % Length of day and night
    
    state = zeros(Nage,Nvacstates,Nvar,4+Nexp,p.Nltlas,length(tvals));
    cases = zeros(Nage,Nvacstates,Nvar,p.Nltlas,length(tvalsi)-1);
    casesI = zeros(Nage,Nvacstates,Nvar,p.Nltlas,length(tvalsi)-1);

    popdat = load('./data/popdata20.mat');
    centdist = sqrt((popdat.centroid(:,1)-popdat.centroid(:,1)').^2 + (popdat.centroid(:,2)-popdat.centroid(:,2)').^2); % Distance between LTLA centroids
    if newvar == 1
        load('./data/vardists.mat');
        impdist = deldist;
    else
        load('./data/impdistdel.mat');
    end
    if surge == 1
        surgestate = zeros(p.Nltlas,tend+delay+27);
    end
    
    if tstart == 279
        l1 = load('./data/Rec1Oct20.mat');
        infatt = 1/160;
        vacdays = getVacDays(tstart,vac,p);
        ICs = getICs(p,Nage,Nvacstates,Nvar,Nexp,vac,vacdays,infatt,l1.Rec,tstart,alpw,alpn,impdistalpha,newvar);
    else
        ICs = prev;
    end
    state(:,:,:,:,:,1) = ICs;
    for i = 1:(length(tvalsi)-1) % Move forward in time
        if ltlaprecmob == 0
            mix = getMix(p, popdat,mob(i),nonworkmov,tvalsi(i),centdist,mobscale,ltlaprecmob); % Get mixing matrix
        else
            mix = getMix(p, popdat,mob(i,:),nonworkmov,tvalsi(i),centdist,mobscale); % Get mixing matrix
        end
        if tstart + i - 1 > deltastart %&& tstart + i - 1 < 508
            tmpplace = randsample(1:312,deln,1,impdist);
            tab = histcounts(tmpplace,0.5:1:312.5);
            pl = find(tab>0);   
            pn = tab(tab>0);
            rem = min(reshape(pn,[1,1,1,1,length(pl)]),state(2,1,1,1,pl,2*(i-1)+1));
            state(2,1,1,1,pl,2*(i-1)+1) = state(2,1,1,1,pl,2*(i-1)+1)-rem;
            state(2,1,Nvar,Nexp+2,pl,2*(i-1)+1) = state(2,1,Nvar,Nexp+2,pl,2*(i-1)+1) + rem;
        end
        for l = 1:2 % Repeat for day and night
            fprec = prec3(:,i,:);
            [foiD,foiN] = getFOI(p,state(:,:,:,:,:,(i-1)*2+l),mix,fprec,agemix,Nage,Nvar,Nexp,Nvacstates); % Force of infection in day and night
            rates = getRates(foiD,foiN,p,state(:,:,:,:,:,(i-1)*2+l),Nage,Nvar,Nexp,Nvacstates,p.Nltlas,l);
            events = getEvents(rates,step(l),state(:,:,:,:,:,(i-1)*2+l),Nage,Nvar,Nexp,Nvacstates,p.Nltlas,neg);
            cases(:,:,:,:,i) = cases(:,:,:,:,i) + squeeze(events(:,:,:,1+Nexp,:)+events(:,:,:,2+Nexp,:));
            casesI(:,:,:,:,i) = casesI(:,:,:,:,i) + squeeze(events(:,:,:,1+Nexp,:));
            state(:,:,:,:,:,(i-1)*2+l+1) = doEventChange(state(:,:,:,:,:,(i-1)*2+l),events,Nexp);
            if l == 1 % Vaccinate in the day
                if surge == 1
                    [state(:,:,:,:,:,(i-1)*2+l+1),surgestate] = doVaccine(state(:,:,:,:,:,(i-1)*2+l+1),p,tvalsi(i),vac,Nexp,Nvar,Nage,Nvacstates,surge,surgeboost,surgethresh,delay,surgestate,omitdd,maxper,svar,threshtype);
                else
                    state(:,:,:,:,:,(i-1)*2+l+1) = doVaccine(state(:,:,:,:,:,(i-1)*2+l+1),p,tvalsi(i),vac,Nexp,Nvar,Nage,Nvacstates,surge);
                end
            end
        end
    end
end
