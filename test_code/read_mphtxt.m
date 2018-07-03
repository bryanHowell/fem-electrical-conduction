function [meshdata] = read_mphtxt(filename)
%% open file

fid=fopen(filename,'r');
if(fid==-1)
    error('The file could not be opened!');
end

%% extract mesh nodes

strtmp='temp';
ck=strfind(strtmp,'number of mesh points');
while( isempty(ck) )
    strtmp=fgets(fid);
    ck=strfind(strtmp,'number of mesh points');
end

meshdata.numnodes=sscanf(strtmp,'%d');
for k=1:3
    fgets(fid);
end

[A,cnt]=fscanf(fid,'%f %f %f',3*meshdata.numnodes);
if(cnt~=meshdata.numnodes*3)
    error('error in reading data');
end

AA=zeros(3,meshdata.numnodes);
AA(:)=A;
AA=AA';
meshdata.nodes=AA;
clear A;
clear AA;

%% extract boundary elements

eleord=2;
if(eleord==2)
    n=3;
    nn_bnd=6;
    nn_int=10;
    sform_bnd='%f %f %f %f %f %f';
    sform_int='%f %f %f %f %f %f %f %f %f %f';
else
    n=1;
    nn_bnd=3;
    nn_int=4;
    sform_bnd='%f %f %f';
    sform_int='%f %f %f %f';    
end

for k=1:n
    strtmp='temp';
    ck=strfind(strtmp,'number of elements');
    while( isempty(ck) )
        strtmp=fgets(fid);
        ck=strfind(strtmp,'number of elements');
    end
end

meshdata.num_bndele=sscanf(strtmp,'%d');
fgets(fid);

[A,cnt]=fscanf(fid,sform_bnd,nn_bnd*meshdata.num_bndele);
if(cnt~=nn_bnd*meshdata.num_bndele)
    error('error in reading data');
end

AA=zeros(nn_bnd,meshdata.num_bndele);
AA(:)=A;
AA=AA';
meshdata.bndele=AA+1; % +1 for 1-n indexing
clear A;
clear AA;

%% extract volume elements

strtmp='temp';
ck=strfind(strtmp,'number of elements');
while( isempty(ck) )
    strtmp=fgets(fid);
    ck=strfind(strtmp,'number of elements');
end

meshdata.num_volele=sscanf(strtmp,'%d');
fgets(fid);

[A,cnt]=fscanf(fid,sform_int,nn_int*meshdata.num_volele);
if(cnt~=nn_int*meshdata.num_volele)
    error('error in reading data');
end

AA=zeros(nn_int,meshdata.num_volele);
AA(:)=A;
AA=AA';
meshdata.volele=AA+1; % +1 for 1-n indexing
clear A;
clear AA;

%% element classification

strtmp='temp';
ck=strfind(strtmp,'number of geometric entity indices');
while( isempty(ck) )
    strtmp=fgets(fid);
    ck=strfind(strtmp,'number of geometric entity indices');
end

fgets(fid);

[A,cnt]=fscanf(fid,'%f',meshdata.num_volele);
if(cnt~=meshdata.num_volele)
    error('error in reading data');
end

meshdata.volele_type=A;

clear A;

%% close file

fclose(fid);