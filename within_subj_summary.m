% within_subj_summary() - Calculate within subject error for single-factor 
% design, using Cousineau (2005) method with Morey (2008) bias-correction.
% 
% Usage:
% >> [NormData,SDev,SErr,CInt] = within_subj_summary(data,[c_val])
% 
% Inputs:
%   data      = 2 or 3-dimensional matrix. Subject (x Time) x Condition.
%               May do multifactorial version in future, see Loftus & Mason
%               (1994).
% Optional inputs:
%   c_val     = Critical value for confidence interval in percent. Default
%               is 95% if not input.
% Outputs:
%   NormData  = Normed data. Subject (x Time) x Condition.
%   SDev      = Standard deviation of within-subject normed data. 
%               (Time) x Condition.
%   SErr      = Standard error of within-subject normed data. df = n-1.
%               (Time) x Condition.
%   CInt      = Confidence interval for the input critical value. Add to
%               and subtract from means to get CI boundaries. 
%               (Time) x Condition.
% 
% Copyright (C) 21/06/2017 Ger Loughnane

function [NormData,SDev,SErr,CInt] = within_subj_summary(varargin)
% read in data
data = varargin{1};
ndim = ndims(data);
% reshape data into 3 dimensional matrix if neccessary, so as to 
% accommodate both 2 and 3-D matrices
if ndim==2
    data = reshape(data,size(data,1),1,size(data,2));
elseif ~ismember(ndim,[2,3])
    error('Data must be 2-D or 3-D matrix')
end
% read in critical value for confidence interval
if nargin==1
    c_val = 0.95;
else
    c_val = varargin{2}/100;
    if c_val<=0 | c_val>=1
        error('Critical value must be > 0% and < 100%')
    end
end

grand_avg = squeeze(mean(mean(data,3),1)); % grand avg
grand_avg = repmat(grand_avg,[size(data,1),1]); % subj
subj_avg = squeeze(mean(data,3)); % subj (optional: x time)

NormData=zeros(size(data)); SDev=zeros(size(data,2),size(data,2)); SErr=zeros(size(data,2),size(data,2)); CInt=zeros(size(data,2),size(data,2));
% iterate through conditions
for i = 1:size(data,3)
    % calculate normed data. new value = old value – subject average + grand average
    NormData(:,:,i) = data(:,:,i)-subj_avg+grand_avg;
    % calculate standard deviation and error from normed data
    SDev(:,i) = std(NormData(:,:,i),[],1);
    SErr(:,i) = squeeze(SDev(:,i)/sqrt(size(data,1)));
    % Confidence interval multiplier for standard error
    % Calculate t-statistic for confidence interval: 
    % e.g., if critical value is .95, use .975 (above/below), and use df=n-1
    tscore = repmat(tinv(c_val/2+0.5,size(data,1)-1),[size(SErr(:,i),1),1]);
    CInt(:,i) = tscore.*SErr(:,i);
end

% Apply correction from Morey (2008) to the standard error and confidence interval
% Get the product of the number of conditions of within-S variables
correction_factor = sqrt(size(data,3)/(size(data,3)-1));
SDev = SDev*correction_factor;
SErr = SErr*correction_factor;
CInt = CInt*correction_factor;

% In case NormData is supposed to be 2-D
NormData = squeeze(NormData);
