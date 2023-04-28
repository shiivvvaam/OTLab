%Simplex
C=[4 10];
A=[2,1;2,5;2,3];
b=[50;100;90];
[m,n]=size(A);
s=eye(m);% identity matrix it creates 
A=[A s b];
cost=zeros(1,n+m+1);
cost(1:n)=C;
bv=n+1:1:n+m;
zjcj=cost(bv)*A-cost;
zcj=[zjcj;A];
simptable=array2table(zcj);
simptable.Properties.VariableNames(1:size(A,2))={'x1','x2','s1','s2','s3','b'};
simptable
flag=true;
% ratio=[];
while flag
  if any(zjcj<0)
        fprintf('curr not optimal');
        zc=zjcj(1:end-1);
        [Enter_val,pvt_col]=min(zc);
    if all(A(:,pvt_col)<=0)
        error('lpp unbounded all entries <=0 in col %d',pvt_col);
    else
        sol=A(:,end);
        col=A(:,pvt_col); 
        for i=1:m
            if col(i)>0
                ratio(i)=sol(i)./col(i);
            else
                ratio(i)=inf;
            end
        end
        [leaving_val,pvt_row]=min(ratio);
    end
    bv(pvt_row)=pvt_col;
    pvt_key=A(pvt_row,pvt_col);
    A(pvt_row,:)=A(pvt_row,:)./pvt_key;
    for i=1:m
        if i~=pvt_row
            A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
        end
    end
zjcj=zjcj-zjcj(pvt_col).*A(pvt_row,:);
zcj=[zjcj;A];
table=array2table(zcj);
table.Properties.VariableNames={'x1','x2','s1','s2','s3','sol'}
else 
        flag=false;
 end
end  
table
