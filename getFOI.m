function [foiD, foiN] = getFOI(p,state,mix,prec,agemix,Nage,Nvar,Nexp,Nvacstates)
%GETFOI Summary of this function goes here
%   Detailed explanation goes here

    % Number of people in day in (age,N)
    dN = reshape(sum(mix,2),[size(mix,1),size(mix,3)]);

    % Number of infected people living in (age, vac, var,N)
    nInf = reshape(state(:,:,:,Nexp+2,:)+p.tau*state(:,:,:,Nexp+3,:),[Nage,Nvacstates,Nvar,p.Nltlas]);

    % Number of infected people in day in (age, vac, var,N)
    dInf = reshape(sum(reshape(mix./reshape(dN,[Nage,1,p.Nltlas]),[Nage,1,1,p.Nltlas,p.Nltlas]).*nInf,4),[Nage,Nvacstates,Nvar,p.Nltlas]);

    foiN = zeros(Nage,Nvacstates,Nvar,p.Nltlas);
    foiD = zeros(Nage,Nvacstates,Nvar,p.Nltlas);
    amix = prec.*agemix;
    for j = 1:Nage
        foiN = foiN + reshape(amix(:,j,:),[size(amix,1),1,1,size(amix,3)]).*(repmat(sum(nInf(j,:,:,:).*p.tt,2),[1,Nvacstates,1,1])./reshape(p.Na(:,j),[1,1,1,p.Nltlas]));
    end
    for j = 1:3
        foiD = foiD + reshape(sum(reshape(amix(:,j,:).*mix./sum(mix,2),[Nage,1,1,p.Nltlas,p.Nltlas]).*(repmat(sum(dInf(j,:,:,:).*p.tt,2),[1,Nvacstates,1,1])./reshape(dN(j,:),[1,1,1,p.Nltlas])),4),[Nage,Nvacstates,Nvar,p.Nltlas]);
    end
    foiN = foiN.*(reshape(p.beta',[1,1,Nvar,p.Nltlas]));
    foiD = foiD.*(reshape(p.beta',[1,1,Nvar,p.Nltlas]));
end


