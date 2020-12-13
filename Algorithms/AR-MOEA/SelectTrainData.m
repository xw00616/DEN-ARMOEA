function [tr_xx, tr_yy] = SelectTrainData(P, N1, N2)

tr_y=P.objs;
NA=size(tr_y,1);
Range = [min(tr_y,[],1);max(tr_y,[],1)];

tr_y = tr_y - repmat(Range(1,:),NA,1);
Cosine = 1 - pdist2(tr_y,tr_y,'cosine');%NA*NA
Cosine(logical(eye(size(Cosine,1)))) = 0;
Choose   = [false(1,NA-N2),true(1,N2)];

while sum(Choose)<N1
        unSelected = find(~Choose);
        [~,x]      = min(max(Cosine(~Choose,Choose),[],2));
        Choose(unSelected(x)) = true;
end

tr_yy=tr_y(Choose,:);
tr_x=P.decs;
tr_xx=tr_x(Choose,:);