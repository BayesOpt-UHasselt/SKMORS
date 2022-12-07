function R = corrMatern(theta,D)
% returns Matern 2.5 Cov Matrix
d = size(theta,1);
k1 = size(D,1);
k2 = size(D,2);

D2=D.^2;

Rfin=ones(k1,k1);

for dd=1:d

Rtemp=(1+ sqrt(5).*D(:,:,dd)./theta(dd)+ 5/3 * D2(:,:,dd)./theta(dd)^2).* exp(- sqrt(5)*D(:,:,dd)./theta(dd));

Rfin=Rfin.*Rtemp;

end

R = Rfin;