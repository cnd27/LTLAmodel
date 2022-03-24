%Set up model script
setUpModel
useFit = 1;
cols = [228,26,28
55,126,184
77,175,74
152,78,163
255,127,0]/255;

if useFit == 1
    % Load fitted parameters
    load('modelPars.mat')
    FF = sum(mean(dists)>0); % Get final generation of fitting
    thetav = parsA(:,:,FF); % Get posterior
else
    thetav = guess;
end
lab= [];
histA = [];
% Get fitted values

    precv = [movmean(thetav(:,precgi(1:36)),4,2) thetav(:,precgi(37:43))];
    mobv = thetav(:,mobgi);
    lab = ["Home working scale","Distance move scale","Sympt. prob. (0-17)" "Sympt. prob. (18-65)" "Sympt. prob. (66+)" "Asympt. trans.","Initial alpha","Daily delta","Rep. (0-17) T1 (%)" "Rep. (0-17) T2 (%)" "Rep. (18-65) (%)" "Rep. (66+) (%)"];
    histA = [ mobv mobscalevdav tauv alpnv delnv repv];
    mobscalev = 3+5*thetav(:,mobscalegi);
    histA = [histA mobscalev];
    dav = thetav(:,transgi(1:3));
    tauv = thetav(:,transgi(4));
    histA = [histA dav tauv];
    alpnv = round(3000*thetav(:,alpngi));
    histA = [histA alpnv];
    delnv = round(20*thetav(:,delngi));
    histA = [histA delnv];
    repv = thetav(:,repgi);
    histA = [histA repv];


figure
for i = 1:length(lab)
    subplot(ceil(length(lab)/ceil(sqrt(length(lab))))+1,ceil(sqrt(length(lab))),i)
    histogram(histA(:,i));
    xlabel(lab{i})
end
subplot(ceil(length(lab)/ceil(sqrt(length(lab))))+1,ceil(sqrt(length(lab))),(ceil(length(lab)/ceil(sqrt(length(lab))))*ceil(sqrt(length(lab)))+1):(ceil(length(lab)/ceil(sqrt(length(lab))))+1)*ceil(sqrt(length(lab))))
precprc = prctile(precv,[2.5 97.5]);
hold on
plot([datetime(2020,1,279) repelem(datetime(2020,1,286:7:573),2) datetime(2020,1,580)],repelem(mean(precv),2),'Color',[cols(2,:)])
pgon1 = fill([datetime(2020,1,279) repelem(datetime(2020,1,286:7:573),2) datetime(2020,1,580) fliplr([datetime(2020,1,279) repelem(datetime(2020,1,286:7:573),2) datetime(2020,1,580)])],...
    [repelem(precprc(1,:),2) fliplr(repelem(precprc(2,:),2))],[cols(2,:)],'LineStyle','none');
set(pgon1,'FaceAlpha',0.4);
ylabel('Reduction in trans.')
ylim([0 1])
xlim([datetime(2020,1,279) datetime(2020,1,580)])
% print(strcat('ParsM',name),'-r300','-dpng')