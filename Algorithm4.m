% use this to test:

% b = [0 10 20 50 80 100 250 500 1000];
% signal = 100*(0.15*exp(-b*0.2) + 0.85*exp(-b*0.001));
% data = signal + 2*randn(size(b));
% [MMSE,MAP,curveFit,logPr,logLh,logPost] = IVIMBayesianEstimation(b,data);
% data = signal + 15*randn(size(b));
% [MMSE,MAP,curveFit,logPr,logLh,logPost] = IVIMBayesianEstimation(b,data);

function Output = Algorithm4(b,data)

% make sure they are column arrays
b = double(b(:));
data = double(data(:));
%class(data)
data1=data(:);
sumData2 = sum(data.^2);

% array of weights in case the NSA at different b-values are different
% weight(i) = NSA for data(i)
weight = ones(size(data));
organ = 'kidney';


weight = weight(:);

% ADC monoexp fit
ln_S=log(data);
p = polyfit(b, ln_S, 1);
ADC = -(p(1));
% range for the parameters
logDgrid = 120; logD = linspace(-10,-4,logDgrid);
logDsgrid = 120; logDs = linspace(-12,0,logDsgrid); % adjusted Ds to be larger range?
fGrid = 60; f = linspace(0,1,fGrid);

% log-priors
if strcmp(organ,'kidney')
    logDmean = -6.2;
    logDstd = 1;
    logPrD = -0.5*(logD - logDmean).^2/logDstd^2 - 0.5*log(2*pi*logDstd^2);
    %
    logDsmean = -3.5;%init: -3.5 
    logDsStd = 1;
    logPrDs = -0.5*(logDs - logDsmean).^2/logDsStd^2 - 0.5*log(2*pi*logDsStd^2);
elseif strcmp(organ,'liver')
    logDmean = -7.0;
    logDstd = 1;
    logPrD = -0.5*(logD - logDmean).^2/logDstd^2 - 0.5*log(2*pi*logDstd^2);
    %
    logDsmean = -3.5;%init: -3.5 
    logDsStd = 1;
    logPrDs = -0.5*(logDs - logDsmean).^2/logDsStd^2 - 0.5*log(2*pi*logDsStd^2);
end
Pr = (exp(logPrD')*exp(logPrDs)).*(repmat(logD',1,logDsgrid)<repmat(logDs,logDgrid,1));
logPr = log(Pr);

% various modelling arrays
exp_b_logD = repmat(exp(-b*exp(logD)),[1 1 logDsgrid]);
exp_b_logDs = permute(repmat(exp(-b*exp(logDs)),[1 1 logDgrid]),[1 3 2]);

% compute log-likelihood and log posterior by looping over f (loops in D
% and D* are incorporated here using array operations)
[logLh,logPost] = deal(zeros(fGrid,logDgrid,logDsgrid));
for nf = 1:fGrid
    G = f(nf)*exp_b_logDs + (1-f(nf))*exp_b_logD;
    logLh(nf,:,:) = squeeze(-0.5*length(b)*log(sumData2 - sum(repmat(data.*weight,[1 logDgrid logDsgrid]).*G).^2./sum(G.^2.*repmat(weight,[1 logDgrid logDsgrid]))));
    logPost(nf,:,:) = logLh(nf,:,:) + reshape(logPr,[1 logDgrid logDsgrid]);
end

% log-posterior
Post = exp(logPost - max(logPost(:)));
% normalisation
sP = sum(Post(:));


% compute MMSE estimates
Output.D = exp(sum(sum(sum(Post,1),3).*logD)/sP);
Output.Ds = exp(sum(squeeze(sum(sum(Post,1),2)).*logDs')/sP);
Output.f = sum(sum(sum(Post,2),3).*f')/sP;


% compute curve for MMSE estimate
G = Output.f*exp(-b*Output.Ds) + (1-Output.f)*exp(-b*Output.D);
curveFit = G*inv(G'*G)*G'*data;
curveFit1 = curveFit/max(curveFit(:));
data=data/max(data(:));

% compute root square error
yresid = data - curveFit1;
SSresid = sum(yresid.^2);
SStotal = (length(data)-1) * var(data);
rsq = 1 - SSresid/SStotal;

Output.Residual=rsq;


end
