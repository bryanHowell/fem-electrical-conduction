function [Ap,fp] = enforce_dirichletbc(A,f)

global meshdata;

%%

ck_dbc=meshdata.bndtype==1;
num_dbc=sum(ck_dbc);

indx_dbc=meshdata.bndnodes(ck_dbc);
val_dbc=meshdata.bndval(ck_dbc);

Ap=A;
fp=f;

%%

for ii=1:num_dbc
    
    jj=indx_dbc(ii);
    fp=fp-Ap(:,jj)*val_dbc(ii);
    
    Ap(jj,:)=0; % clear jjth column
    Ap(:,jj)=0; % clear jjth row
    Ap(jj,jj)=1; % set LHS as a known
    
    fp(jj)=val_dbc(ii);
    
end