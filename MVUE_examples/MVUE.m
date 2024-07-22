%% MVUE as described by Cohn et al., 1989, WRR

% Parameters:
% p - the regression coefficients of the rating curve output by polyval, i.e. [p,~] = polyval(log(Q), log(Qs),1)
% W - Daily mean discharge (Q) values
% Q - sampled discharge values (i.e., same as those used to calculate p)
% Qs - sampled Qs values (i.e., same as those used to calculate p)
% center - the center of the data after Cohn et al. (1992)
% Returns:
% yld - the mean yeilds for each value of W (presumably daily or hourly averaged discharge)
% ci_h - the upper 95% CI bound of yld 
% ci_l - the lower 95% CI bound of yld


% Tested against loadest with Taiwan sed data, LOADEST AMLE mean yields should yield exactly the
% same results here, except in Taiwan assume an input of cms instead of cfs, and
% output of metric instead of short tons (why does USGS use cfs and short ton? god only knows).
% The implementation of the CIs here is based on Gilroy et al., 1990 and
% differs from the implementation in LOADEST - there are 
% differences in the CIs from loadest of up to 0.3 percent sometimes, and I cannot pin down why.


function [yld,ci_l,ci_h]=MVUE(p,Q,Qs,W,center)
    %Initialize
    n1=length(Q);
    n2=length(W);
    ld=log(Q); 
    X=[ones(n1,1) ld(:)]; 
    lw=log(W); 
    Xb=[ones(n2,1) lw(:)]; 
    % Calculate the covariance matrix Vt between W and Q
    Vt=Xb/(X'*X)*Xb';
    % Extract the diagonals - variances
    V=zeros(1,length(Xb));
    for i=1:length(Vt)
        V(i)=Vt(i,i);
    end
    % Sum of squared residuals 
    SE1=sum((log(Qs)-(p(2)+(log(Q)-center).*p(1))).^2)./(length(Qs)-2);
    %Degrees of freedom
    m=length(Q)-2;
    
    % MLE factor
    fact=(m+1)/(2*m)*((1-V)*SE1);
    
    % Initial value for the series expansion
    bt=fact.*m^2/(2*(m+1));

    gm1=1;
    term1=1;
    gmi=gm1;
    % Perform iterative series expansion to calculate the MLE
    for r=1:1000
        term1=term1.*bt./((m/2+r-1).*r);
        gm1=gm1+term1;
        if abs(gmi-gm1)<1e-10
            break
        end
    end
    % Calculate the estimated yield for each
    b1= polyval([p(1), p(2)],log(W)-center);
    b1= exp(b1);    
    % the final yield
    yld=b1(:).*gm1(:);
    
    %
    % Calculate variance
    var=0;
    % Iterate over asset pairs 
    for i=1:length(lw(:))
            for j=1:i
                % Calculate the covariance term
                fact=(((1-V(i))*(1-V(j))/4*SE1.^2));
    
                bl=fact;
                gm2=1;
                term1=1;
                gmi=gm2;
    
                % Iteration from Cohn et al., 1989
                for r=1:1000
                    term1=term1.*bl./((m/2+r-1).*r);
                    gm2=gm2+term1;
                    if abs(gmi-gm2)<1e-20
                        break
                    end
                    gmi=gm2;
                end
                % Calculate the Term (Gilroy et al., 1990) for the variance
                Term = (((b1(i)*b1(j)*exp((V(i)+V(j)+2*Vt(i,j)) *SE1/2)) * ...
                    (exp((2-V(i)-V(j))*SE1/2)) * ...
                    (gm2)) - ...
                    (exp(log(b1(i))+(SE1/2)) * ...
                    exp(log(b1(j))+(SE1/2))));
                % Update the variance 
                if i==j
                    var=var+Term;
                else
                    var=var+2*Term;
                    
                end
            end
    end
    for i = 1:length(yld)
        var = var + yld(i).^2 * (exp(SE1)-1); 
    end
    %Standard error of the mean
    U=sqrt(var)/length(lw(:));

    B=sqrt(log(1+(U/mean(yld))^2));

    % log mean
    A=log(mean(yld))-B^2/2;
    % Upper and lower bounds of the confidence interval
    ci_h=exp(A+1.96*B);
    ci_l=exp(A-1.96*B);
end
