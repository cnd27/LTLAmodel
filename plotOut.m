%% Model set up
p = getParameters(); % Get model parameters
Nvar = 3; % Number of variants
Nvacstates = 3; %Number of vaccine states (unvac, v1, v2)
Nexp = 1; % Number of compartments in exposed class
Nage = size(p.Na,2); % Number of age classes
agemixg = getAgeMix(); % Get age mixing matrix

load('./data/NewVac14.mat') % Load vaccine data (delayed by 14 days)
vac = vac14;
load('./data/googleWork.mat') % Load Google working data
nonworkmov = 0; % 0 if non-workers do not move between LTLAs

load('./data/LineListData.mat') % Load linelist data (with and without age breakdown)
load('./data/LineListDataAge.mat')
dataW =  permute(casesLTLAwa(40:82,:,:),[3,2,1]); % Consider weekly data 5 October 2020 - 1 August 2021
dataW(dataW==0) = 1; %Remove zeros to avoid NaNs

load('./data/SgeneLTLA.mat') % Load S gene data
sgtf = num(:,35:335)./tot(:,35:335); % Consider daily data 5 October 2020 - 1 August 2021
sgtfa = sum(num(:,35:335))./sum(tot(:,35:335)); % Data summed over LTLAs
% Replace NaNs with local average
av = 0;
while sum(isnan(sgtf),'all') > 0
    newmov = movmean(num(:,35:335),av+1,2)./movmean(tot(:,35:335),av+1,2);
    sgtf(isnan(sgtf)) = newmov(isnan(sgtf));
    av = av+1;
end

% Precautionary behaviour and mobility derived from Google data
load('./data/googleMov.mat')
gprec = (100+movav(tstart:(tend-1),:))/100;
gmob = max(0,min(1,1-(100+movwork(tstart:(tend-1),:))/100));

neg = 0; % Negative binmoial for exposures
surge = 0; % Do surge vaccine
surgeboost = 0; % Boost vaccine by percentage
surgethresh = 0; % LTLA threshold for surge vaccination
impdistalpha = 1;
newvar = 1;
ltlaprecmob = 1;
load('./data/regionnum.mat')
printing = 0;
plotmaps = 0;

W = 43; % Length of data in weeks
tstart = 279; % Start in days since 1st Jan 2020
tend = 580; % End in days

% Number of fitted parameters
Npars = 56;
precgi = 1:43;
mobgi = 44;
mobscalegi = 45;
transgi = 46:49;
alpngi = 50;
delngi = 51;
repgi = 52:55;

% Load fitted parameters
load('modelPars.mat')
FF = sum(mean(dists)>0); % Get final generation of fitting
thetav = parsA(:,:,FF); % Get posterior

% Get fitted values

precv = [movmean(thetav(:,precgi(1:36)),4,2) thetav(:,precgi(37:43))];
mobv = thetav(:,mobgi);
mobscalev = 3+5*thetav(:,mobscalegi);
dav = thetav(:,transgi(1:3));
tauv = thetav(:,transgi(4));
R0v = 2.7.*ones(size(parsA,1),1);
trv = [4/2.7 6.5/4].*ones(size(parsA,1),1);
prec3v = repmat(precv,[1,1,3]);
agemixv = repmat(reshape(agemixg,[1,3,3]),[size(parsA,1),1,1]);
alpnv = round(3000*thetav(:,alpngi)); 
alpwv = 2*ones(size(parsA,1),1);
delnv = round(20*thetav(:,delngi));
deltastartv = 426*ones(size(parsA,1),1);
reglv = ones(size(parsA,1),312);   

%% Model Simulation

% Run simulations
NT = 8; % Number of realisations
prev = 1; % Alternative ICs (ignore)
% Get cases 
cases = zeros(Nage,Nvacstates,Nvar,p.Nltlas,tend-tstart,NT);
casesI = zeros(Nage,Nvacstates,Nvar,p.Nltlas,tend-tstart,NT);
casesrep = zeros(Nage,Nvacstates,Nvar,p.Nltlas,tend-tstart,NT);
casesIrep = zeros(Nage,Nvacstates,Nvar,p.Nltlas,tend-tstart,NT);
states = zeros(Nage,Nvacstates,Nvar,4+Nexp,p.Nltlas,1+2*(tend-tstart));
parfor i = 1:NT
    fprintf('%d\n',i)
    idx = randi(size(thetav,1)); % Choose random posterior set
    prec = repelem(precv(idx,:),7)';
    prec3 = squeeze(prec3v(idx,:,:));
    prec3 = reshape(repelem(prec3(:),7),[tend-tstart,Nage])';
    mob = mobv(idx)*gmob;
    mobscale = mobscalev(idx,:);
    da = dav(idx,:);
    tau = tauv(idx);
    R0 = R0v(idx);
    tr = trv(idx,:);
    agemix = squeeze(agemixv(idx,:,:));
    alpn = alpnv(idx);
    alpw = alpwv(idx);
    deln = delnv(idx);
    deltastart = deltastartv(idx);
    regl = reglv(idx,:);
    [state,casest,casesIt] = runLTLAmodel(p,vac,Nvar,Nexp,Nvacstates,Nage,nonworkmov,tstart,tend,prev,prec3,mob,agemix,da,tau,R0,tr,alpw,alpn,deln,mobscale,impdistalpha,deltastart,neg,newvar,ltlaprecmob,surge,surgeboost,surgethresh); %[279 580]
    states = states + state;
    repo = repmat(reshape([repelem(thetav(idx,repgi(1)),22*7), repelem(thetav(idx,repgi(2)),21*7);repelem(thetav(idx,repgi(3)),43*7);repelem(thetav(idx,repgi(4)),43*7)],[Nage,1,1,1,tend-tstart]),[1,Nvacstates,Nvar,p.Nltlas,1]);
    casesIrep(:,:,:,:,:,i) = binornd(casesIt,repo);
    casesrep(:,:,:,:,:,i) = casesIrep(:,:,:,:,:,i) + binornd(casest-casesIt,repo);
end
states = states/NT;

meancasesrep = mean(squeeze(sum(sum(sum(casesrep,2),3),4)),3);
meancasesrepw = squeeze(sum(reshape(meancasesrep,[3,7,43]),2));
meanSGTF = mean(squeeze(sum(sum(sum(casesrep(:,:,2,:,:,:))),4))./squeeze(sum(sum(sum(sum(casesrep))))),2);

prccasesrep = prctile(squeeze(sum(sum(sum(sum(casesrep)))))',[2.5,97.5]);
prccasesrepage = zeros(2,301,3);
for i = 1:3
    prccasesrepage(:,:,i) = prctile(squeeze(sum(sum(sum(casesrep(i,:,:,:,:,:)))))',[2.5,97.5]);
end
prcsgtf = prctile((squeeze(sum(sum(sum(casesrep(:,:,2,:,:,:))),4))./squeeze(sum(sum(sum(sum(casesrep))))))',[2.5,97.5]);
prccasesrepagew = zeros(2,43,3);
for i = 1:3
    prccasesrepagew(:,:,i) = prctile(squeeze(sum(reshape(squeeze(sum(sum(sum(casesrep(i,:,:,:,:,:))))),[7,43,NT])))',[2.5,97.5]);
end

cols = [228,26,28
55,126,184
77,175,74
152,78,163
255,127,0]/255;

%% Plot model output

h = figure();
set(h, 'PaperUnits','centimeters')
set(h, 'PaperSize',[24 8])
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
subplot(1,3,1)
hold on
a = {};
plot(datetime(2020,1,tstart:(tend-1)),casesd(tstart:(tend-1)),'.','Color',[cols(1,:)])
a = plot(datetime(2020,1,tstart:(tend-1)),sum(meancasesrep),'Color',[cols(1,:)]);
pgon1 = fill([datetime(2020,1,tstart:(tend-1)) fliplr(datetime(2020,1,tstart:(tend-1)))],[prccasesrep(1,:) fliplr(prccasesrep(2,:))],[cols(1,:)],'LineStyle','none');
set(pgon1,'FaceAlpha',0.4);
set(gca,'FontSize',8,'FontName','Arial')
ylabel('Daily cases','FontSize',8,'FontName','Arial')

subplot(1,3,2)
hold on
plot(datetime(2020,1,tstart:(tend-1)),sgtfa,'.','Color',[cols(1,:)])
plot(datetime(2020,1,tstart:(tend-1)),meanSGTF,'Color',[cols(1,:)])
pgon2 = fill([datetime(2020,1,tstart:(tend-1)) fliplr(datetime(2020,1,tstart:(tend-1)))],[prcsgtf(1,:) fliplr(prcsgtf(2,:))],[cols(1,:)],'LineStyle','none');
set(pgon2,'FaceAlpha',0.4);
set(gca,'FontSize',8,'FontName','Arial')
ylabel('Proportion of cases with SGTF','FontSize',8,'FontName','Arial')

subplot(1,3,3)
hold on
a = {};
for i = 1:3
    plot(datetime(2020,1,tstart:(tend-1)),casesda(tstart:(tend-1),i),'.','Color',[cols(i,:)])
    a{i} = plot(datetime(2020,1,tstart:(tend-1)),meancasesrep(i,:),'Color',[cols(i,:)]);
    pgon1 = fill([datetime(2020,1,tstart:(tend-1)) fliplr(datetime(2020,1,tstart:(tend-1)))],[prccasesrepage(1,:,i) fliplr(prccasesrepage(2,:,i))],[cols(i,:)],'LineStyle','none');
    set(pgon1,'FaceAlpha',0.4);
end
leg = legend([a{1} a{2} a{3}],{'0-17','18-65','66+'});
set(gca,'FontSize',8,'FontName','Arial')
set(leg,'FontSize',8,'FontName','Arial')
ylabel('Daily cases','FontSize',8,'FontName','Arial')
if printing == 1
    print(strcat('DailyM',name),'-r300','-dpng')
end

load('./data/popdata20.mat')
h = figure();
set(h, 'PaperUnits','centimeters')
set(h, 'PaperSize',[36 17])
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
chooseLTLAs = zeros(7,3);
compL = sum(sum((permute(sum(reshape(squeeze(sum(sum(mean(casesrep,6),2),3)),[3,312,7,43]),3),[4,2,1,3])-casesLTLAwa(40:82,:,:)).^2,1),3);
for l = 1:7
    
    nn = find(regionnum==l);
    [~,li] = sort(compL(nn));
    orderLTLA = nn(li);
    chooseLTLAs(l,1) = orderLTLA(1);
    chooseLTLAs(l,2) = orderLTLA(floor(length(nn)/2));
    chooseLTLAs(l,3) = orderLTLA(end);
    
    for ll = 1:3
        subplot(3,7,l+(ll-1)*7)
        hold on
        a = {};
        mmm = mean(squeeze(sum(sum(casesrep(:,:,:,chooseLTLAs(l,ll),:,:),2),3)),3);
        prcmmm = zeros(2,301,3);
        for i = 1:3
            prcmmm(:,:,i) = prctile(squeeze(sum(sum(casesrep(i,:,:,chooseLTLAs(l,ll),:,:))))',[2.5,97.5]);
        end
        for i = 1:3
            title(string(ltlaname(chooseLTLAs(l,ll))),'FontSize',8,'FontName','Arial');
            plot(datetime(2020,1,tstart:(tend-1)),casesLTLAa(tstart:(tend-1),chooseLTLAs(l,ll),i),'.','Color',[cols(i,:)])
            a{i} = plot(datetime(2020,1,tstart:(tend-1)),mmm(i,:),'Color',[cols(i,:)]);
            pgon1 = fill([datetime(2020,1,tstart:(tend-1)) fliplr(datetime(2020,1,tstart:(tend-1)))],[prcmmm(1,:,i) fliplr(prcmmm(2,:,i))],[cols(i,:)],'LineStyle','none');
            set(pgon1,'FaceAlpha',0.4);
        end
        set(gca,'FontSize',8,'FontName','Arial')
        if l == 7 && ll == 1
            leg = legend([a{1} a{2} a{3}],{'0-17','18-65','66+'});
            set(leg,'FontSize',8,'FontName','Arial')
            ylabel('Daily cases','FontSize',8,'FontName','Arial')
        end
    end
end
if printing == 1
    print(strcat('LTLAspreadM',name),'-r300','-dpng')
end

h = figure();
set(h, 'PaperUnits','centimeters')
set(h, 'PaperSize',[16 8])
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
subplot(1,2,1)
hold on
a = {};
for i = 1:3
    a{i} = plot([datesw(40) repelem(datesw(41:82),2) datesw(83)],repelem(sum(casesLTLAwa(40:82,:,i),2),2),'Color',[cols(i,:)]);
    plot([datesw(40) repelem(datesw(41:82),2) datesw(83)],repelem(meancasesrepw(i,:),2),':','Color',[cols(i,:)]);
    pgon1 = fill([datesw(40) repelem(datesw(41:82),2) datesw(83) fliplr([datesw(40) repelem(datesw(41:82),2) datesw(83)])],[repelem(prccasesrepagew(1,:,i),2) fliplr(repelem(prccasesrepagew(2,:,i),2))],[cols(i,:)],'LineStyle','none');
    set(pgon1,'FaceAlpha',0.4);
end
leg = legend([a{1} a{2} a{3}],{'0-17','18-65','66+'});
set(gca,'FontSize',8,'FontName','Arial')
set(leg,'FontSize',8,'FontName','Arial')
ylabel('Weekly cases','FontSize',8,'FontName','Arial')
subplot(1,2,2)
hold on
a = {};
for i = 1:3
    plot([datesw(41) repelem(datesw(42:82),2) datesw(83)],repelem(sum(casesLTLAwa(41:82,:,i),2)./sum(casesLTLAwa(40:81,:,i),2),2),'Color',[cols(i,:)])
    a{i} = plot([datesw(41) repelem(datesw(42:82),2) datesw(83)],repelem(meancasesrepw(i,2:43)./meancasesrepw(i,1:42),2),':','Color',[cols(i,:)]);

end
set(gca,'FontSize',8,'FontName','Arial')
ylabel('Change in weekly cases','FontSize',8,'FontName','Arial')
if printing == 1
    print(strcat('WeeklyM',name),'-r300','-dpng')
end
%% Plot maps

if plotmaps == 1
    pltla = shaperead('/Users/christopherdavis/OneDrive - University of Warwick/JUNIPER/Spatial/LTLA2020/plotLTLAsimple.shp');
    keytimes = [279 306 336 367 398 426 457 487 518 548 573];
    keytimesW = [keytimes;keytimes+1;keytimes+2;keytimes+3;keytimes+4;keytimes+5;keytimes+6];
    dataCase = permute(sum(reshape(casesLTLA(keytimesW(:),:),[7,length(keytimes),p.Nltlas])),[3,2,1]);
    dataCase(dataCase==0)=1;
    modelCase = squeeze(sum(squeeze(sum(sum(mean(casesrep(:,:,:,:,(keytimes-tstart+1),:),6),2),3))));
    modelCaseDist = squeeze(sum(squeeze(sum(sum(casesrep(:,:,:,:,(keytimes-tstart+1),:),2),3))));
    for q = 1:6
        modelCase = modelCase + squeeze(sum(squeeze(sum(sum(mean(casesrep(:,:,:,:,(keytimes-tstart+1+q),:),6),2),3))));
        modelCaseDist = modelCaseDist + squeeze(sum(squeeze(sum(sum(casesrep(:,:,:,:,(keytimes-tstart+1+q),:),2),3))));
    end
    
    modelDiff = modelCase./dataCase;
    modelDiff(modelDiff==0)=1;
    percentile = 100*abs(NT/2-sum(modelCaseDist<dataCase,3))/(NT/2); 
    percentile2 = 100*abs(NT/2-sum(modelCaseDist./sum(modelCaseDist)<dataCase./sum(dataCase),3))/(NT/2);
    h = figure();
    set(h, 'PaperUnits','centimeters')
    set(h, 'PaperSize',[80 22])
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
        
    for k = 1:length(keytimes)
        subplot(4,11,k)
        hold on
        title(datestr(datetime(2020,1,keytimes(k))))
        mapDynScaled = 1+round(250*max(0,min(4,log10(dataCase(:,k)))));
        cmap = flipud(spring(1001));
    
        axis off
        for i = 1:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        ylim([49 57])
        colormap(gca,cmap)
        if k == 11
            cc = colorbar;
            cc.Label.String = 'Data cases';
            cc.Location = 'eastoutside';
            cc.Label.FontName = 'Arial';
            cc.Label.FontSize = 8;
            cc.FontSize = 8;
            cc.Ticks = 0:0.25:1;
            cc.TickLabels = [1 10 100 1000 10000];
            cc.Position = cc.Position + [0.05 -0.05 0 0];
        end
    
        caxis([0 1])
        AxesHandle=findobj(gcf,'Type','axes');
            a1 = get(AxesHandle(1), 'Position');
    
        axes('units','normalized','position',[a1(1)-0.0085 a1(2)-0.01 0.06*0.6*0.65 0.25*0.4*0.65],'Color','none')
        londontext2 = text(-0.49,51.775,'London','VerticalALignment','top','FontSize',6,'FontName','Arial');
        box on
        set(gca,'XTick',[],'YTick',[])
        ylim([51.286 51.78])
        hold on
        for i = 281:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        subplot(4,11,k+11)
        hold on
        mapDynScaled = 1+round(250*max(0,min(4,log10(modelCase(:,k)))));
        cmap = flipud(spring(1001));
    
        axis off
        for i = 1:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        ylim([49 57])
        colormap(gca,cmap)
        if k == 11
            cc = colorbar;
            cc.Location = 'eastoutside';
            cc.Label.String = 'Mean model cases';
            cc.Label.FontName = 'Arial';
            cc.Label.FontSize = 8;
            cc.FontSize = 8;
            cc.Ticks = 0:0.25:1;
            cc.TickLabels = [1 10 100 1000 10000];
            cc.Position = cc.Position + [0.05 -0.05 0 0];
        end
    
        caxis([0 1])
        AxesHandle=findobj(gcf,'Type','axes');
            a1 = get(AxesHandle(1), 'Position');
    
        axes('units','normalized','position',[a1(1)-0.0085 a1(2)-0.01 0.06*0.6*0.65 0.25*0.4*0.65],'Color','none')
        londontext2 = text(-0.49,51.775,'London','VerticalALignment','top','FontSize',6,'FontName','Arial');
        box on
        set(gca,'XTick',[],'YTick',[])
        ylim([51.286 51.78])
        hold on
        for i = 281:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        subplot(4,11,k+22)
        hold on
        mapDynScaled = 501+round(500*max(-1,min(1,log10(modelDiff(:,k)))));
        cmap = flipud(parula(1601));
        cmap = cmap([1:430 731:871 1172:1601],:);
    
        axis off
        for i = 1:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        ylim([49 57])
        colormap(gca,cmap)
        if k == 11
            cc = colorbar;
            cc.Label.String = 'Model difference';
            cc.Location = 'eastoutside';
            cc.Label.FontName = 'Arial';
            cc.Label.FontSize = 8;
            cc.FontSize = 8;
            cc.Ticks = 0.5+0.5*log10([0.1 0.2 0.5 1 2 5 10]);
            cc.TickLabels = [0.1 0.2 0.5 1 2 5 10];
            cc.Position = cc.Position + [0.05 -0.05 0 0];
        end
    
        caxis([0 1])
        AxesHandle=findobj(gcf,'Type','axes');
            a1 = get(AxesHandle(1), 'Position');
    
        axes('units','normalized','position',[a1(1)-0.0085 a1(2)-0.01 0.06*0.6*0.65 0.25*0.4*0.65],'Color','none')
        londontext2 = text(-0.49,51.775,'London','VerticalALignment','top','FontSize',6,'FontName','Arial');
        box on
        set(gca,'XTick',[],'YTick',[])
        ylim([51.286 51.78])
        hold on
        for i = 281:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        subplot(4,11,k+33)
        hold on
        mapDynScaled = 1+round(50*(max(80,percentile2(:,k))-80));
        cmap = flipud(bone(501));
        cmap =[repmat(cmap(1,:),[500,1]);cmap];
    
        axis off
        for i = 1:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        ylim([49 57])
        colormap(gca,cmap)
        if k == 11
            cc = colorbar;
            cc.Label.String = 'Percent interval that contains data';
            cc.Location = 'eastoutside';
            cc.Label.FontName = 'Arial';
            cc.Label.FontSize = 8;
            cc.FontSize = 8;
            cc.Ticks = 0:0.25:1;
            cc.TickLabels = 80:5:100;
            cc.Position = cc.Position + [0.05 -0.05 0 0];
        end
    
        caxis([0 1])
        AxesHandle=findobj(gcf,'Type','axes');
            a1 = get(AxesHandle(1), 'Position');
    
        axes('units','normalized','position',[a1(1)-0.0085 a1(2)-0.01 0.06*0.6*0.65 0.25*0.4*0.65],'Color','none')
        londontext2 = text(-0.49,51.775,'London','VerticalALignment','top','FontSize',6,'FontName','Arial');
        box on
        set(gca,'XTick',[],'YTick',[])
        ylim([51.286 51.78])
        hold on
        for i = 281:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
    end
    pause(2)
    AxesHandle=findobj(gcf,'Type','axes');
    ap = get(AxesHandle, 'Position');
    
    for i = 2:2:88
        set(AxesHandle(i), 'Position',ap{i}+[-0.018 -0.12 0.02 0.13]);
    end
    if printing == 1
        print(strcat('MapDataCase2M',name),'-r300','-dpng')
    end
    
    dataCaseA = permute(sum(reshape(casesLTLAa(keytimesW(:),:,:),[7,length(keytimes),p.Nltlas,3])),[3,2,4,1]);
    dataCaseA(dataCaseA==0)=1;
    modelCaseA = permute(squeeze(sum(sum(mean(casesrep(:,:,:,:,(keytimes-tstart+1),:),6),2),3)),[2,3,1]);
    modelCaseDistA = permute(squeeze(sum(sum(casesrep(:,:,:,:,(keytimes-tstart+1),:),2),3)),[2,3,1,4]);
    for q = 1:6
        modelCaseA = modelCaseA + permute(squeeze(sum(sum(mean(casesrep(:,:,:,:,(keytimes-tstart+1+q),:),6),2),3)),[2,3,1]);
        modelCaseDistA = modelCaseDistA + permute(squeeze(sum(sum(casesrep(:,:,:,:,(keytimes-tstart+1+q),:),2),3)),[2,3,1,4]);
    end
    
    modelDiffA = modelCaseA./dataCaseA;
    modelDiffA(modelDiffA==0)=1;
    percentileA = 100*abs(NT/2-sum(modelCaseDistA<dataCaseA,4))/(NT/2); 
    h = figure();
    set(h, 'PaperUnits','centimeters')
    set(h, 'PaperSize',[60 22])
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
        
    for k = 1:length(keytimes)
        subplot(3,11,k)
        hold on
        title(datestr(datetime(2020,1,keytimes(k))))
        mapDynScaled = 1+round(50*(max(80,percentileA(:,k,1))-80));
        cmap = flipud(bone(501));
        cmap =[repmat(cmap(1,:),[500,1]);cmap];
    
        axis off
        for i = 1:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        ylim([49 57])
        colormap(gca,cmap)
        if k == 11
            cc = colorbar;
            cc.Label.String = 'Percent interval that contains data';
            cc.Location = 'eastoutside';
            cc.Label.FontName = 'Arial';
            cc.Label.FontSize = 8;
            cc.FontSize = 8;
            cc.Ticks = 0:0.25:1;
            cc.TickLabels = 80:5:100;
            cc.Position = cc.Position + [0.05 -0.05 0 0];
        end
    
        caxis([0 1])
        AxesHandle=findobj(gcf,'Type','axes');
            a1 = get(AxesHandle(1), 'Position');
    
        axes('units','normalized','position',[a1(1)-0.0085 a1(2)-0.01 0.06*0.6*0.65 0.25*0.4*0.65],'Color','none')
        londontext2 = text(-0.49,51.775,'London','VerticalALignment','top','FontSize',6,'FontName','Arial');
        box on
        set(gca,'XTick',[],'YTick',[])
        ylim([51.286 51.78])
        hold on
        for i = 281:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        subplot(3,11,k+11)
        hold on
        mapDynScaled = 1+round(50*(max(80,percentileA(:,k,2))-80));
        cmap = flipud(bone(501));
        cmap =[repmat(cmap(1,:),[500,1]);cmap];
    
        axis off
        for i = 1:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        ylim([49 57])
        colormap(gca,cmap)
        if k == 11
            cc = colorbar;
            cc.Label.String = 'Percent interval that contains data';
            cc.Location = 'eastoutside';
            cc.Label.FontName = 'Arial';
            cc.Label.FontSize = 8;
            cc.FontSize = 8;
            cc.Ticks = 0:0.25:1;
            cc.TickLabels = 80:5:100;
            cc.Position = cc.Position + [0.05 -0.05 0 0];
        end
    
        caxis([0 1])
        AxesHandle=findobj(gcf,'Type','axes');
            a1 = get(AxesHandle(1), 'Position');
    
        axes('units','normalized','position',[a1(1)-0.0085 a1(2)-0.01 0.06*0.6*0.65 0.25*0.4*0.65],'Color','none')
        londontext2 = text(-0.49,51.775,'London','VerticalALignment','top','FontSize',6,'FontName','Arial');
        box on
        set(gca,'XTick',[],'YTick',[])
        ylim([51.286 51.78])
        hold on
        for i = 281:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        subplot(3,11,k+22)
        hold on
        mapDynScaled = 1+round(50*(max(80,percentileA(:,k,3))-80));
        cmap = flipud(bone(501));
        cmap =[repmat(cmap(1,:),[500,1]);cmap];
    
        axis off
        for i = 1:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
        ylim([49 57])
        colormap(gca,cmap)
        if k == 11
            cc = colorbar;
            cc.Label.String = 'Percent interval that contains data';
            cc.Location = 'eastoutside';
            cc.Label.FontName = 'Arial';
            cc.Label.FontSize = 8;
            cc.FontSize = 8;
            cc.Ticks = 0:0.25:1;
            cc.TickLabels = 80:5:100;
            cc.Position = cc.Position + [0.05 -0.05 0 0];
        end
    
        caxis([0 1])
        AxesHandle=findobj(gcf,'Type','axes');
            a1 = get(AxesHandle(1), 'Position');
    
        axes('units','normalized','position',[a1(1)-0.0085 a1(2)-0.01 0.06*0.6*0.65 0.25*0.4*0.65],'Color','none')
        londontext2 = text(-0.49,51.775,'London','VerticalALignment','top','FontSize',6,'FontName','Arial');
        box on
        set(gca,'XTick',[],'YTick',[])
        ylim([51.286 51.78])
        hold on
        for i = 281:312
            mapshow(pltla(i),'FaceColor',cmap(mapDynScaled(i),:))
        end
    end
    pause(2)
    AxesHandle=findobj(gcf,'Type','axes');
    ap = get(AxesHandle, 'Position');
    
    for i = 2:2:66
        set(AxesHandle(i), 'Position',ap{i}+[-0.018 -0.12 0.02 0.13]);
    end
    if printing == 1
        print(strcat('MapAgeM',name),'-r300','-dpng')
    end
end