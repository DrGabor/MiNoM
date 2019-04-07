function [ y ] = eppdf( data, sita, p )
lamda = p*sita^(1/p)/(2*gamma(1/p));
y = lamda * exp(-sita*(abs(data)).^p);
end

