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
    load(strcat('fitpars',name,'NVP6.mat'))
    FF = sum(mean(dists)>0); % Get final generation of fitting
    thetav = parsA(:,:,FF); % Get posterior
else
    thetav = guess;
end
lab= [];
histA = [];
% Get fitted values
if fitprec == 1
    if smoothprec == 1
%         precv = movmean(thetav(:,precgi),4,2); % Smooth over 4 weeks
        precv = [movmean(thetav(:,precgi(1:36)),4,2) thetav(:,precgi(37:43))];
    else
        precv = thetav(:,precgi);
    end
else 
    precv = thetav(:,precgi);
end
if fitmob == 1
    mobv = thetav(:,mobgi);
    lab = [lab "Home working scale"];
    histA = [histA mobv];
else
    mobv = ones(size(parsA,1),1);
end
if fitmobscale == 1
    mobscalev = 3+5*thetav(:,mobscalegi);
    lab = [lab "Distance move scale"];
    histA = [histA mobscalev];
else
    mobscalev = (1e6)*ones(size(parsA,1),3);
end
if fittrans == 1
    dav = thetav(:,transgi(1:3));
    tauv = thetav(:,transgi(4));
    R0v = 2.7.*ones(size(parsA,1),1);
    trv = [4/2.7 6.5/4].*ones(size(parsA,1),1);
    prec3v = repmat(precv,[1,1,3]);
    lab = [lab "Sympt. prob. (0-17)" "Sympt. prob. (18-65)" "Sympt. prob. (66+)" "Asympt. trans."];
    histA = [histA dav tauv];
elseif fittrans == 2
    dav = thetav(:,transgi(1:3));
    tauv = thetav(:,transgi(4));
    R0v = 2.7.*ones(size(parsA,1),1);
    trv = [4/2.7 6.5/4].*ones(size(parsA,1),1);
    prec3v = reshape([thetav(:,transgi(5)) thetav(:,transgi(6)) ones(size(parsA,1),1)],[size(parsA,1),1,Nage]).*repmat(precv,[1,1,3]);
    lab = [lab "Sympt. prob. (0-17)" "Sympt. prob. (18-65)" "Sympt. prob. (66+)" "Asympt. trans." "Suscept. (0-17)"  "Suscept. (18-65)"];
    histA = [histA dav tauv thetav(:,transgi(5)) thetav(:,transgi(6))];
elseif fittrans == 3
    dav = thetav(:,transgi(1:3));
    tauv = thetav(:,transgi(4));
    R0v = 2 + thetav(:,transgi(7));
    trv = 1.5+0.5*thetav(:,transgi(8:9));
    prec3v = reshape([thetav(:,transgi(5)) thetav(:,transgi(6)) ones(size(parsA,1),1)],[size(parsA,1),1,Nage]).*repmat(precv,[1,1,3]);
    lab = [lab "Sympt. prob. (0-17)" "Sympt. prob. (18-65)" "Sympt. prob. (66+)" "Asympt. trans." "Suscept. (0-17)"  "Suscept. (18-65)" "R0" "Alpha adv." "Delta adv."];
    histA = [histA dav tauv thetav(:,transgi(5)) thetav(:,transgi(6)) R0v trv];
else
    dav = [0.05 0.2 0.7].*ones(size(parsA,1),1);
    tauv = 0.25.*ones(size(parsA,1),1);
    R0v = 2.7.*ones(size(parsA,1),1);
    trv = [4/2.7 6.5/4].*ones(size(parsA,1),1);
    prec3v = repmat(precv,[1,1,3]);
end
if fitagemix == 1
    agemixv = 15*reshape(thetav(:,agemixgi),[size(parsA,1),3,3]);
    lab = [lab "A(1,1)" "A(2,1)" "A(3,1)" "A(1,2)" "A(2,2)" "A(3,2)" "A(1,3)" "A(2,3)" "A(3,3)"];
    histA = [histA reshape(agemixv,[100,9])];
else
    agemixv = repmat(reshape(agemixg,[1,3,3]),[size(parsA,1),1,1]);
end
if fitalpn == 1
    alpnv = round(3000*thetav(:,alpngi));
    lab = [lab "Initial alpha"];
    histA = [histA alpnv];
else
    alpnv = 3000*fitalpng*ones(size(parsA,1),1);
end
if fitalpw == 1
    alpwv = thetav(:,alpwgi);
    lab = [lab "Alpha weighting"];
    histA = [histA alpwv];
else
    alpwv = 5*fitalpwg*ones(size(parsA,1),1);
end
if fitdeln == 1
    delnv = round(20*thetav(:,delngi));
    lab = [lab "Daily delta"];
    histA = [histA delnv];
else
    delnv = 20*fitdelng*ones(size(parsA,1),1);
end
if fitdeltastart == 1
    deltastartv = 415+round(50*thetav(:,deltastartgi));
    lab = [lab "Delta start"];
    histA = [histA deltastartv];
else
    deltastartv = deltastart*ones(size(parsA,1),1);
end
if fitreg == 1
    regv = thetav(:,reggi);
    reglv = zeros(size(parsA,1),312);
    for z = 1:7
        reglv(:,regionnum==z) = repmat(regv(:,z),[1,sum(regionnum==z)]);
    end
    lab = [lab "London" "East of England" "South East" "South West" "Midlands" "North West" "North East"];
    histA = [histA regv];
    regionpop = zeros(1,7);
    for z = 1:7
        regionpop(z) = sum(p.N(regionnum==z))./sum(p.N);
    end
    precv = precv.*sum(regv.*regionpop,2);
else
    reglv = ones(size(parsA,1),312);
end
if fitrep == 1
    repv = thetav(:,repgi);
    lab = [lab "Rep. (0-17) T1 (%)" "Rep. (0-17) T2 (%)" "Rep. (18-65) (%)" "Rep. (66+) (%)"];
    histA = [histA repv];
else
    repv = ones(size(parsA,1),4);
end

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
print(strcat('ParsM',name),'-r300','-dpng')