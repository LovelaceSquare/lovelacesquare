function [r2,lofc]=lofNaN(de,c,s)
dr=c*s;
res=de-dr;
de=de(:); res=res(:);
iN=find(~isnan(res));
lofc=100*sqrt((res(iN)'*res(iN))/(de(iN)'*de(iN)));
r2=(1-(lofc/100)^2)*100;
end