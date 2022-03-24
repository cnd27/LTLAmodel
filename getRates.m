function rates = getRates(foiD,foiN,p,state,nage,nvar,nexp,nvacstate,NN,l)
%GETRATES Rate for each event

rates = zeros(nage,nvacstate,nvar,4+nexp,NN);

if l == 1
    foi = foiD;
else
    foi = foiN;
end
sigmaT = repmat(p.sigma,[1,1,nvar,1,NN]);
dT = repmat(p.d,[1,1,nvar,1,NN]);
rates(:,:,:,1,:) = reshape(foi,[nage,nvacstate,nvar,1,p.Nltlas]).*sigmaT.*repmat(sum(state(:,:,:,1,:),3),[1,1,nvar,1,1]);
rates(:,:,:,1+nexp,:) = dT.*p.epsilon.*state(:,:,:,nexp+1,:);
rates(:,:,:,2+nexp,:) = (1-dT).*p.epsilon.*state(:,:,:,nexp+1,:);
rates(:,:,:,3+nexp,:) = p.gamma.*state(:,:,:,nexp+2,:);
rates(:,:,:,4+nexp,:) = p.gamma.*state(:,:,:,nexp+3,:);

if sum(rates<0,'all') > 0
    error('Something went wrong with rates')
end

end

