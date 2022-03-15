load('./data/numtot2.mat')
load('./data/totf.mat')
load('./data/popdata20.mat')
NHS=readtable('/Users/christopherdavis/OneDrive - University of Warwick/JUNIPER/UKCoV/data/locs.csv');

region = cell(312,1);
for i = 1:length(ltlaname)
    region{i} = NHS{find(ltlaname(i)==NHS{:,2}),1};
end
region{65} = 'North West';
region{73} = 'Midlands';
region{171} = 'East of England';
region{176} = 'North West';
region{193} = 'South East';
R = strings(312,1);
for i = 1:312
    R(i) = string(region{i});
end

rSW = R == 'South West';
rSE = R == 'South East';
rE = R == 'East of England';
rL = R == 'London';
rNW = R == 'North West';
rNE = R == 'North East and Yorkshire';
rM = R == 'Midlands';
r = [rSW rSE rE rL rNW rNE rM];

NHScases = zeros(7,371);
NHScasesS = zeros(7,371);
NHScasesSF = zeros(7,371);
for i = 1:7
    NHScases(i,:) = sum(totf(r(:,i),:));
    NHScasesS(i,:) = sum(tot(r(:,i),:));
    NHScasesSF(i,:) = sum(num(r(:,i),:));
end

alphaT = [268 275   282   289   296   303   310   317   324   331   338   345   352   359   366   373   380   387   394   401   408   415   422   429]-244;
deltaT = [419 426   433   440   447   454   461   468   475   482   489   496   503   510   517   524   531   538   545   552   559   566   573   580   587   594]-244; 

AlphaCases = zeros(7,length(alphaT)-2);
AlphaCasesChange = zeros(7,length(alphaT)-2);
AlphaPos = zeros(7,length(alphaT)-2);
AlphaPosTest = zeros(7,length(alphaT)-2);
for i = 1:size(AlphaCases,2)
    AlphaCases(:,i) = sum(NHScases(:,alphaT(i+1):(alphaT(i+2)-1)),2);
    AlphaCasesChange(:,i) = sum(NHScases(:,alphaT(i+1):(alphaT(i+2)-1)),2)./sum(NHScases(:,alphaT(i):(alphaT(i+1)-1)),2);
    AlphaPos(:,i) = sum(NHScasesSF(:,alphaT(i+1):(alphaT(i+2)-1)),2);
    AlphaPosTest(:,i) = sum(NHScasesS(:,alphaT(i+1):(alphaT(i+2)-1)),2);
end
AlphaSGTF = AlphaPos./AlphaPosTest;

DeltaCases = zeros(7,length(deltaT)-2);
DeltaCasesChange = zeros(7,length(alphaT)-2);
DeltaPos = zeros(7,length(deltaT)-2);
DeltaPosTest = zeros(7,length(deltaT)-2);
for i = 1:size(DeltaCases,2)
    DeltaCases(:,i) = sum(NHScases(:,deltaT(i+1):(deltaT(i+2)-1)),2);
    DeltaCasesChange(:,i) = sum(NHScases(:,deltaT(i+1):(deltaT(i+2)-1)),2)./sum(NHScases(:,deltaT(i):(deltaT(i+1)-1)),2);
    DeltaPos(:,i) = sum(NHScasesSF(:,deltaT(i+1):(deltaT(i+2)-1)),2);
    DeltaPosTest(:,i) = sum(NHScasesS(:,deltaT(i+1):(deltaT(i+2)-1)),2);
end
DeltaSGTF = DeltaPos./DeltaPosTest;