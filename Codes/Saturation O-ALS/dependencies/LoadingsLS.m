function [P] = LoadingsLS(D,T)
[~,J]=size(D);
P=nan(size(T,2),size(D,2));
for j=1:J
    ireal=find(isfinite(D(:,j)));
    if isempty(ireal)
        P(:,j)=nan;
    else
        P(:,j)=T(ireal,:)\D(ireal,j);
    end
end
end

