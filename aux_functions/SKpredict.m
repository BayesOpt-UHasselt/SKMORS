function [f, DT_MSE]  = SKpredict(model,Xpred,Bpred)

% edited by Sebastian 20/06/2018:
%    * MSE is also calculated: 
%       - DT_MSE= deterministic MSE

% make predictions at prediction points using a stochastic kriging model  
% model = output of SKfit
% Xpred = (K x d) matrix of prediction points
% Bpred = (K x b) matrix of basis functions at each prediction point
%         The first column must be a column of ones!
% f = (K x 1) predictions at predictions points
% 
% Exmaples
%      SK_gau  = SKpredict(skriging_model,XK,ones(K,1));
% Based on parameter estimates skriging_model obtained from SKfit.m,
% use SK model to predict the values at prediction points XK with constant
% prediction trend, Bpred = ones(K,1)

% retrieve model parameters from model structure obtained from SKfit
X = model.X;
minX = model.minX;
maxX = model.maxX;
[k d] = size(X);
theta = model.theta;
gammaP = model.gamma;
beta = model.beta;
Z2 = model.Z2;
% L = model.L;
tau2 = model.tausquared;
sigmahatinv= model.sigmahatinv;
sigmaDTinv= model.sigmahatinvDT;

% simple check for dimensions of Xpred and X
K = size(Xpred,1);     % number of prediction points
if (size(Xpred,2)~=d)
    error('Prediction points and design points must have the same dimension (number of columns).');
end
if (size(Bpred,1)~=K)
    error('Basis function and prediction point matrices must have the same number of rows.')
end
if not(all(Bpred(:,1)==1))
    error('The first column of the basis function matrix must be ones.')
end

% calculate distance matrix for prediction points
Xpred = (Xpred - repmat(minX,K,1)) ./ repmat(maxX-minX,K,1);
if gammaP == 2
    distXpred =  abs(repmat(reshape(Xpred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K])).^2;
else
    distXpred =  abs(repmat(reshape(Xpred', [1 d K]),[k,1,1]) ...
        - repmat(X,[1 1 K]));
end

% calculate correlations between prediction points and design points
D = distXpred;
if gammaP == 3
    T = repmat(reshape(theta,[1 d 1]),[k 1 K]);
    Sigm_pred = tau2*prod(((D<=(T./2)).*(1-6*(D./T).^2+6*(D./T).^3) ...
        +((T./2)<D & D<=T).*(2*(1-D./T).^3)),2);
elseif gammaP == 4
    Sigm_pred = tau2.*corrMatern_pred(theta,D);
else
    Sigm_pred = tau2*exp(sum(-D.*repmat(reshape(theta,[1 d]),[k 1 K]),2));
end

Sigm_pred = reshape(Sigm_pred,[k K 1]);

f= Bpred*beta + Sigm_pred'*sigmahatinv*Z2; % stochastic kriging prediction

DT_MSE=zeros(K,1);

for ci=1:K
 delta= 1-ones(k,1)'* sigmaDTinv * Sigm_pred(:,ci);
 DT_MSE(ci)= tau2-Sigm_pred(:,ci)' * sigmaDTinv * Sigm_pred(:,ci) + delta^2/(ones(k,1)'* sigmaDTinv * ones(k,1));
end


