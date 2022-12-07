function R = corrMatern_pred(theta,D)
% return matern correlation function
d = size(theta,1);

D2=D.^2;

Rfin=ones(size(D,1),1,size(D,3));

for dd=1:d

Rtemp=(1+ sqrt(5).*D(:,dd,:)./theta(dd)+ 5/3 .* D2(:,dd,:)./theta(dd)^2).* exp(- sqrt(5)*D(:,dd,:)./theta(dd));

Rfin=Rfin.*Rtemp;

end

R = Rfin;