function [Precision, Validation_Results_Array] = drugCIPHER_MS_Validation(Drug_Target_Relation, PPI_Adjacent_Matrix, Drug_Similarity_Matrix1, Square1, Drug_Similarity_Matrix2, Square2, Repeat_Times, Validation_Set_Size)


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
%%if nargin < 5
%%    Repeat_Times = 100;
%%    Validation_Set_Size = 20;
%%elseif nargin == 5
%%    Validation_Set_Size = 20;
%%end


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


%% Leave One Out Validation
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

Drug_Similarity_Matrix1_data = Drug_Similarity_Matrix1_data.^Square1;

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

Drug_Similarity_Matrix2_data = Drug_Similarity_Matrix2_data.^Square2;

disp('Computing Leave One Out Validation...');
Precision(Repeat_Times) = 0;
Validation_Results_Array(Repeat_Times, Drug_Target_Relation_Number) = 0;
ValidationSet(Validation_Set_Size) = 0;
ValidationSet_Gene_score(Validation_Set_Size) = 0;

for TT = 1:Repeat_Times
    fprintf('    The %dth Interation in Validation Procedure...\n',TT);
    success = 0;
    TotalTimes = 0;
                    ValidationSet(Validation_Set_Size) = 0;
                ValidationSet_Gene_score(Validation_Set_Size) = 0;
    for i = 1:DrugNum %For all drugs, TargetsSet is the validation set; the number i-th drug

        RealTargets = Drug2Targets{i};
        
        %fprintf('\t\t The %dth Drug Candidates...\n',i);
        if isempty(RealTargets)
            continue; %next drug;
        else
            [nunt, RealTargetsetSize] = size(RealTargets);
            RealTargets = RealTargets + 1;
            for rr = 1:RealTargetsetSize
                
                Gene2Drug_Closeness_Tempt = Gene2Drug_Closeness;
                
                ValidationSet(Validation_Set_Size) = 0;
                ValidationSet_Gene_score(Validation_Set_Size) = 0;
                
                %r=unidrnd(k);  % randomly generate one known target for the drug;
                SelectedTarget = RealTargets(rr); % RealTarget is the real target in the validation set
                ValidationSet(1) = SelectedTarget; % ValidationSet: the target set with one real target in it;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%subtract the
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%removed gene
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%from
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Target2Drug
                
                Vector = exp(-PPI_Shortest_Dist(SelectedTarget,:).^2);
                Gene2Drug_Closeness_Tempt(:,i) = Gene2Drug_Closeness(:,i)-Vector';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ii = 2;
                %while ii<=RealTargetsetSize
                while ii <= Validation_Set_Size
                    r = unidrnd(GeneNum);
                    Redundency = 0;  %examing whether the new randomly generated targets is in the realtargets set;
                    for jj = 1:RealTargetsetSize
                        if RealTargets(jj) == r
                            Redundency = 1;
                            break;
                        end
                    end
                    if Redundency ~= 1
                        ValidationSet(ii) = r;
                        ii=ii+1;
                    else
                        continue;
                    end
                end
                %%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%sort the target in ValidationSet and find the first target whether
                %%%%%%%%%%%%%%%%the real target;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%20%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                Drug_i_similarity2D = Drug_Similarity_Matrix1_data(i,:)';
                Drug_i_similarityATC = Drug_Similarity_Matrix2_data(i,:)';
                Similarity(:,1) = ones(size(Drug_i_similarity2D));
                Similarity(:,2) = Drug_i_similarity2D;
                var2 = var(Drug_i_similarity2D);
                Similarity(:,3) = Drug_i_similarityATC;
                var3 = var(Drug_i_similarityATC);
                
                TotalTimes = TotalTimes+1;
                
                Gene_In_Validation_Set2Drugs = Gene2Drug_Closeness_Tempt(ValidationSet,:)';
                
                Rho2 = corr(Gene_In_Validation_Set2Drugs,Similarity(:,2));
                Rho3 = corr(Gene_In_Validation_Set2Drugs,Similarity(:,3));
                Beta = Similarity\Gene_In_Validation_Set2Drugs;
                ValidationSet_Gene_score = (sqrt(var2)./abs(Beta(3,:)).*Rho2' + sqrt(var3)./abs(Beta(2,:)).*Rho3')./sqrt(var2./Beta(3,:).^2 + var3./Beta(2,:).^2);
                
                if ValidationSet_Gene_score(1) == max(ValidationSet_Gene_score)
                    success = success + 1;
                    Validation_Results_Array(TT,TotalTimes) = 1;
                else
                    Validation_Results_Array(TT,TotalTimes) = 0;
                end
                %%%%%%%%%%%%%%%%
            end
        end
    end
    Precision(TT) = success/TotalTimes;
end
