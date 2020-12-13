function xnew=IndividualSelect(PopDec, PopObj, PopMSE, Ke, flag)

    [IDX, ~]=kmeans(PopObj,Ke);
    UIDX=unique(IDX);  
    
    if flag==0
        Uncertainty=mean(PopMSE,2);
        for i=1:length(UIDX)
            Pindex=find(IDX==UIDX(i));
            [~,ind]=max(Uncertainty(Pindex));
            %disp(sprintf('ind is %u', ind)); disp(Pindex');
            xnew(i,:)=PopDec(Pindex(ind),:);          
        end
    else
        Convergence=sqrt(sum(PopObj.^2,2));
        for i=1:length(UIDX)
            Pindex=find(IDX==UIDX(i));
            [~,ind]=min(Convergence(Pindex));
            %disp(sprintf('ind is %u', ind));disp(Pindex');
            xnew(i,:)=PopDec(Pindex(ind),:);      
        end        
    end

    
