%The historical "center" of Q as described by Cohn et al., 1992
%Q - historic Q values
% returns c, the (logged) geometric center of the data

function c= calc_center(Q)
    mn=mean(log(Q(Q>0)));
    Q=(log(Q(Q>0)));
    c =  mn + sum((Q-mn).^3)/sum((2*(Q-mn).^2));
end