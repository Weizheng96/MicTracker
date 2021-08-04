function [mu, sigma] = fitTruncGauss(data)

% data = clsstIntDiff;
% data = clsstOvrlpRate;
% figure; histogram(data,100); xlim([min(data),max(data)]);


% muDif = mean(data);
muDif = median(data);
sigmaDif = std(data);
% figure; hist(normrnd(muDif,sigmaDif,[10^5,1]),200); xlim([min(data),max(data)]);

if (true)
    mu0 = muDif;
    sigma0 = sigmaDif;
    % pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma) ./ (normcdf((mu0+sigma0),mu,sigma)-normcdf((mu0-sigma0),mu,sigma));
    start = [mu0,sigma0];
    truncedVec = data(data>=(mu0-sigma0) & data<=(mu0+sigma0));
    timesCount = 1;
    while var(truncedVec) == 0
        timesCount  = timesCount + 1;
        truncedVec = data(data>=(mu0-timesCount*sigma0) & data<=(mu0+timesCount*sigma0));
    end
    pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma) ./ (normcdf((mu0+timesCount*sigma0),mu,sigma)-normcdf((mu0-timesCount*sigma0),mu,sigma));
    % figure; hist(truncedVec,10); xlim([min(data),max(data)]);
    [paramEsts,~] = mle(truncedVec, 'pdf',pdf_truncnorm, 'start',start, 'lower',[-Inf 0]);
    mu = paramEsts(1);
    sigma = paramEsts(2);
    % % figure; hist(normrnd(mu,sigma,[10^5,1]),200); xlim([min(data),max(data)]);
    
else
    mu = muDif;
    sigma = sigmaDif;
end

% %%%  Display
% figure; hd = histogram(data,200,'Normalization','probability'); xlim([min(data),max(data)]);
% tmpX = (hd.BinEdges(2:201) + hd.BinEdges(1:200)) / 2;
% tmp = hist(normrnd(mu,sigma,[10^5,1]),tmpX);
% tmp = tmp/sum(tmp);
% hold on; plot(tmpX,tmp, 'LineWidth', 1.5);

