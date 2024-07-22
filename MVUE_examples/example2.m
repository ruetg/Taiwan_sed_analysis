load('example_data');
Q=d;
C=c;
%C is measured sediment concentration in units of ppm, Q is measured water discharge m3/s, W contains daily
%averaged water discharge m3/s

% In many cases, there are fewer daily discharge
% measurements than anticipated, or the daily discharge is 0 kg.  In Taiwan, rarely does this mean that the discharge
% is actually 0, insteamd this is more likely a problem with the readings.
% In this case the data is not censored in the traditional sense, and we do not account for censored data here.

I=(Q>0&C>0);
Q = Q(I);
C = C(I);
W = W(W>0);


% Before doing the fit I converted measured sediment discharge at each station to kg/m^3.  
% (multiply by 86.4 (kg s) / (m3 d ppm)). 
% If your data is in different format from m3/s, you will have to adjust this coefficient.  
% MVUE function is agnostic to timescale so it doesn't matter if your discharge 
% is integrated over daily or hourly averaged Q values, 
% just make sure this coefficient is reflective of the format your data.
coef_ppm = 86.4;

Qs =Q .* C * coef_ppm;

% Fit Q-Qs data using polyfit

% We regress over the log-centered data. When comparing a values between years, for example, 
% the "center" should be based on Q from all the available years, not
% simply the year of interest. Thus you would have to modify the "Q" value
% that goes into the center_calc operation here if you are selecting from a
% larger superset of data and comparing values of a between them. For
% example the "center" of all the data for this station (S7) in Taiwan is
% actually 4.1296, higher than here. % This centering should only affect the regression coefficient a, 
% and not the yields.
Q_center = calc_center(Q);
[p,~]=polyfit(log(Q) - Q_center,log(Qs),1);

%compute the discharge based on W 
[yld1,ci_l,ci_u]=MVUE(p, Q, Qs, W, Q_center);

%Visualize the output and "centered" Q
loglog(W,yld1,'.')
xlabel('Q (daily mean) m^3/s')
ylabel('Q_s (kg / day)')
hold on;
plot([exp(Q_center),exp(Q_center)+1e-20],[min(yld1),max(yld1)])
[~,S]=polyfit(log(W),log(yld1),1);

%This determines the quality of the fit - Does MVUE yeild a nearly log-linear fit? 
%We cannot use MVUE if the Q-Qs relationship is not sufficiently strong, or if there is too much
%extraplation with W beyond the range of Q. For this we adopted the ad-hoc
%approach of finding the r2 value between Q and predicted Qs, to make sure
%it is above a threshold that ensures it is mostly a log-linear
%relationship.
r2 = 1 - (S.normr/norm(log(yld1) - mean(log(yld1))))^2;

if r2<.99
    disp('Warning: This data does not yeild a good fit and should not be used with MVUE')
end

% The total yield will be the sum of the daily-averaged Qs.  
% We extrapolate based on the mean discharge for the given time period, i.e. 1 month.


Total_annual_yield = mean(yld1) * 365 / 1e9; %kg to megatons andp 365.25 days per year
upper_ci = ci_u * 365.25 / 1e9;
lower_ci = ci_l * 365.25 / 1e9;

disp(['Total annual yield is: ', num2str(Total_annual_yield), ' Mt', ' with 95% CI ', ...
  num2str(lower_ci), '-',num2str(upper_ci),' mT'])


%%% Rating parameter b is generally reported as the slope of log(C) vs log(Q),
%  whereas we did the regression on Qs, so we must subtract 1 from the exponent
rating_b = p(1) - 1;

% Rating parameter ã is affected by log-transformation bias, so we cannot directly report the regression coefficient
%  Instead we report the concentration at the centered discharge value calculated using MVUE:


[yld_center, ~, ~]=MVUE(p,Q,Qs,[exp(Q_center),1],Q_center);

rating_a = yld_center(1) / exp(Q_center) / coef_ppm; %MVUE yield at center, divided by Q at center, gives concentration

disp(['The corresponding rating parameters are ǎ = ', num2str(rating_a) , ' ppm and b = ', num2str(rating_b)])











