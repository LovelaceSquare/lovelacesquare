function [C] = ScoresLS(D,S)
[I,~]=size(D);
C=nan(size(D,1),size(S,1));
for i=1:I
    jreal=find(isfinite(D(i,:)));
    if isempty(jreal)
        C(i,:)=nan;
    else
        C(i,:)=D(i,jreal)/S(:,jreal);
    end
end
end
