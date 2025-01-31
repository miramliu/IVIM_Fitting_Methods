% use this to test:

% b = [0 10 20 50 80 100 250 500 1000];
% signal = 100*(0.15*exp(-b*0.2) + 0.85*exp(-b*0.001));
% data = signal + 2*randn(size(b));
% [MMSE,MAP,curveFit,logPr,logLh,logPost] = IVIMBayesianEstimation(b,data);
% data = signal + 15*randn(size(b));
% [MMSE,MAP,curveFit,logPr,logLh,logPost] = IVIMBayesianEstimation(b,data);

function MMSE = IVIMBayesianEstimation_simulation(b,data,organ,patientID,display,weight)

% make sure they are column arrays
b = double(b(:));
data = double(data(:));
%class(data)
data1=data(:);
sumData2 = sum(data.^2);

% array of weights in case the NSA at different b-values are different
% weight(i) = NSA for data(i)
if nargin<=5
    weight = ones(size(data));
end
if nargin<=4
    %display = 1;
    display = 0;
end
if nargin<=3
    patientID = 0;
end
if nargin<=2
    %organ = 'liver';
    %organ = 'brain'; %%??
    organ = 'kidney';
end

weight = weight(:);

% ADC monoexp fit
ln_S=log(data);
p = polyfit(b, ln_S, 1);
ADC = -(p(1));
% range for the parameters
logDgrid = 120; logD = linspace(-10,-4,logDgrid);
logDsgrid = 120; logDs = linspace(-100,0,logDsgrid); % adjusted Ds to be larger range?
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
MMSE.D = exp(sum(sum(sum(Post,1),3).*logD)/sP);
%MMSE.Dva = exp(sum(sum(sum(Post,1),3).*(logD-MMSE.D).^2)/sP);
MMSE.Ds = exp(sum(squeeze(sum(sum(Post,1),2)).*logDs')/sP);
MMSE.f = sum(sum(sum(Post,2),3).*f')/sP;
IVIM_data=[MMSE.D,MMSE.f, MMSE.Ds, ADC];
% find MAP estimate (this will be close the the LS estimate)
[~,I] = max(max(max(Post,[],3)));
MAP.D = exp(logD(I));
[~,I] = max(max(max(Post,[],2),[],1),[],3);
MAP.Ds = exp(logDs(I));
[~,I] = max(max(max(Post,[],2),[],3));
MAP.f = f(I);

% compute curve for MMSE estimate
G = MMSE.f*exp(-b*MMSE.Ds) + (1-MMSE.f)*exp(-b*MMSE.D);
curveFit = G*inv(G'*G)*G'*data;
curveFit1 = curveFit/max(curveFit(:));
data=data/max(data(:));

% compute root square error
yresid = data - curveFit1;
SSresid = sum(yresid.^2);
SStotal = (length(data)-1) * var(data);
rsq = 1 - SSresid/SStotal;

if display==1
%if 1 == 1; %hardcoding to show
    % display the curve and the marginal posterior distribution
    g1=figure(1);
    pos = get(g1,'position');
%     set(gcf,'name',sprintf('Patient id #%d',patientID));
    set(g1,'position',[pos(1:2) 605  806])
    title(char(patientID));
    subplot(3,1,1); plot(b,data,'o',b,curveFit1,'r'); xlabel('b'); ylabel('Signal')
    subplot(3,2,3); plot(logD,sum(sum(Post,1),3)); xlabel('log D'); ylabel('pdf')
    subplot(3,2,4); plot(logDs,squeeze(sum(sum(Post,1),2))); xlabel('log D*'); ylabel('pdf')
    subplot(3,2,5); plot(f,sum(sum(Post,2),3)); xlabel('f'); ylabel('pdf')
    subplot(3,2,6); contour(logDs,logD,squeeze(sum(Post,1)),100); xlabel('log D*'); ylabel('log D'); title('joint pdf')
    %savefig(g1,['/Users/octaviabane/Dropbox/MATLAB/RenalTx/RenalTX_IVIM/' patientID '.fig'], 'compact');
%     savefig(['/Users/octaviabane/Dropbox/MATLAB/Diffusion/RenalTX_IVIM/' char(patientID) '.fig'], 'compact');
     
                  
    g2=figure(2); 
    plot(b, data1, 'o', b, curveFit, 'r'); xlabel('bval'); ylabel('Signal');
    %savefig(g2,['/Users/octaviabane/Dropbox/MATLAB/RenalTx/RenalTX_IVIM/' patientID '_full_signal.fig'], 'compact');

%     savefig([char(patientID) '_full_signal.fig']);
    pause();

end    
