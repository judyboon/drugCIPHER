function [Drug_Gene_Score_Matrix] = drugCIPHER_MS_Overall(Drug_Target_Relation, PPI_Adjacent_Matrix, Drug_Similarity_Matrix1, square1, Drug_Similarity_Matrix2, square2)





%%%%%%%%%%

%% Drug_Target_Relation:
%%% Drug_Targe_Relation contains targets information for each drug, every
%%% line represents a drug, and the corresponding known targets are
%%% seperated by '\t'. If the line is '-1', there is no known target for
%%% that drug. The targets are represented by the gene index in the PPI
%%% network. The index are begin from 0. 


%% PPI_Adjacent_Matrix
%%% PPI_Adjacent_Matrix is the adjacent matrix for the input PPI network.
%%% 1 presents having edge; otherwise is 0

%% Drug_Similarity_Matrix1, Drug_Similarity_Matrix2
%%% The drug similarity matrix, corresponding to the index in Drug_Target_Relation.

%% Repeat_Times
%%% The random repeat times of the validation procedure, the default value
%%% is 100;



%% Set the default value


%% Begin Alogrithm
%%%%% 
%%%%%
%%%%%

%% Compute Shortest Distance in the PPI Network
disp('Computing The Shortest Distance Between Nodes...');


Suffix_Exame = regexp(PPI_Adjacent_Matrix,'\w*txt','match');
tag = 0;
if isempty(Suffix_Exame)
    tag = 1;
    Suffix_Exame = regexp(PPI_Adjacent_Matrix,'\w*mat','match');
    if isempty(Suffix_Exame)
        disp('error, please check the input format');
        return;
    end
end

if tag == 0
    PPI_Adjacent_Matrix_data = load(PPI_Adjacent_Matrix);
else
    Struct_Data = load(PPI_Adjacent_Matrix);
    Field_Name = fieldnames(Struct_Data);
    PPI_Adjacent_Matrix_data = getfield(Struct_Data,Field_Name{1});
end
%PPI_Adjacent_Matrix_data(find(PPI_Adjacent_Matrix_data == -1)) = 0;
Sparse_Matrix = sparse(PPI_Adjacent_Matrix_data);
clear PPI_Adjacent_Matrix_data;
%Sparse_Matrix = Sparse_Matrix + Sparse_Matrix';
%Sparse_Matrix(Sparse_Matrix > 1) = 1;
PPI_Shortest_Dist = graphallshortestpaths(Sparse_Matrix,'Directed',false);
%PPI_Shortest_Dist(isnan(PPI_Shortest_Dist)) = inf;
[GeneNum,t] = size(PPI_Shortest_Dist);

%% Compute Gene2drug Closeness
disp('Computing Gene to Drug Closeness Matrix...');

Suffix_Exame = regexp(Drug_Target_Relation,'\w*txt','match');
if isempty(Suffix_Exame)
    disp('error, please check the input format');
    return;
end

fid = fopen(Drug_Target_Relation);
Drug_Target_Relation_Number = 0;
DrugNum = 1;
line = fgetl(fid);
while ischar(line)
    tArray = regexp(line,'\t','split');
    if ~isempty(line)
        Array = [];
        [t,TargetNum] = size(tArray);
        if strcmp(tArray(1),'-1')
            Drug2Targets{DrugNum} = [];
        else
            Drug_Target_Relation_Number = Drug_Target_Relation_Number + TargetNum;
            Array(TargetNum) = 0;
            for i = 1:TargetNum
                Array(i) = str2num(tArray{i});
            end           
            Drug2Targets{DrugNum} = Array;            
        end
        line = fgetl(fid);
        DrugNum = DrugNum + 1;
    else
        break;
    end
end
DrugNum = DrugNum - 1;
fclose(fid);

Gene2Drug_Closeness(GeneNum,DrugNum) = 0;

for i = 1:GeneNum
    for j = 1:DrugNum
        Array = Drug2Targets{j};
        if ~isempty(Array)
            Gene2Drug_Closeness(i,j) = sum(exp(-(PPI_Shortest_Dist(i,Array+1)).^2));
        else
            Gene2Drug_Closeness(i,j) = 0;
        end
    end
end

%%
disp('Loading Drug Similarity Matrixes...');

%%%%
%%%% check input format for Similarity_Matrix1
Suffix_Exame = regexp(Drug_Similarity_Matrix1,'\w*txt','match');
tag = 0;
if isempty(Suffix_Exame)
    tag = 1;
    Suffix_Exame = regexp(Drug_Similarity_Matrix1,'\w*mat','match');
    if isempty(Suffix_Exame)
        disp('error, please check the input format');
        return;
    end
end

if tag == 0
    Drug_Similarity_Matrix1_data = load(Drug_Similarity_Matrix1);
else
    Struct_Data = load(Drug_Similarity_Matrix1);
    Field_Name = fieldnames(Struct_Data);
    Drug_Similarity_Matrix1_data = getfield(Struct_Data,Field_Name{1});
end

Drug_Similarity_Matrix1_data = Drug_Similarity_Matrix1_data.^square1;

%%%%
%%%% check input format for Similarity_Matrix2
Suffix_Exame = regexp(Drug_Similarity_Matrix2,'\w*txt','match');
tag = 0;
if isempty(Suffix_Exame)
    tag = 1;
    Suffix_Exame = regexp(Drug_Similarity_Matrix2,'\w*mat','match');
    if isempty(Suffix_Exame)
        disp('error, please check the input format');
        return;
    end
end

if tag == 0
    Drug_Similarity_Matrix2_data = load(Drug_Similarity_Matrix2);
else
    Struct_Data = load(Drug_Similarity_Matrix2);
    Field_Name = fieldnames(Struct_Data);
    Drug_Similarity_Matrix2_data = getfield(Struct_Data,Field_Name{1});
end

Drug_Similarity_Matrix2_data = Drug_Similarity_Matrix2_data.^square2;

Drug_Gene_Score_Matrix(DrugNum,GeneNum) = 0;

for i = 1:DrugNum %For all drugs, TargetsSet is the validation set; the number i-th drug
    fprintf('    Computing Target Profile for the %dth drug...\n',i);
    Drug_i_similarity2D = Drug_Similarity_Matrix1_data(i,:)';
    Drug_i_similarityATC = Drug_Similarity_Matrix2_data(i,:)';
    Similarity(:,1) = ones(size(Drug_i_similarity2D));
    Similarity(:,2) = Drug_i_similarity2D;
    var2 = var(Drug_i_similarity2D);
    Similarity(:,3) = Drug_i_similarityATC;
    var3 = var(Drug_i_similarityATC);
    %Mcov=sum((Similarity(:,2)-mean(Similarity(:,2))).*(Similarity(:,3)-mean(Similarity(:,3))))/DrugNum;
    %Mcov=(Similarity(:,2),Similarity(:,3));
    %for ii=1:RealTargetsetSize
    
    Rho2 = corr(Gene2Drug_Closeness',Similarity(:,2));
    Rho3 = corr(Gene2Drug_Closeness',Similarity(:,3));
    Beta = Similarity\Gene2Drug_Closeness';
    Drug_Gene_Score_Matrix(i,:) = (sqrt(var2)./abs(Beta(3,:)).*Rho2' + sqrt(var3)./abs(Beta(2,:)).*Rho3')./sqrt(var2./Beta(3,:).^2 + var3./Beta(2,:).^2);
    
    %%%%%%%%%%%%%%%%
end


