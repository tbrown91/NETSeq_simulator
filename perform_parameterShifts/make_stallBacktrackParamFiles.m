load ../Scer/NETseq_geneNames.mat;
load ../compare_paramsAllGenes/NETseqBestFitBinned.mat;
load ../compare_paramsAllGenes/NETseqGeneParamsBinned.mat;

model_idx = 7;

%% Make files of parameters to test
double_stallFile = fopen('WT_NETseqDoubleStallParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==3)||(i==7)
            fprintf(double_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*2);
        else
            fprintf(double_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(double_stallFile,'\n');
end
fclose(double_stallFile);

%% 
triple_stallFile = fopen('WT_NETseqTripleStallParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==3)||(i==7)
            fprintf(triple_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*3);
        else
            fprintf(triple_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(triple_stallFile,'\n');
end
fclose(triple_stallFile);

%% 
quadruple_stallFile = fopen('WT_NETseqQuadrupleStallParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==3)||(i==7)
            fprintf(quadruple_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*4);
        else
            fprintf(quadruple_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quadruple_stallFile,'\n');
end
fclose(quadruple_stallFile);

%% Make files of parameters to test
double_backtrackFile = fopen('WT_NETseqDoubleBacktrackParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==5)||(i==9)
            fprintf(double_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*2);
        else
            fprintf(double_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(double_backtrackFile,'\n');
end
fclose(double_backtrackFile);

%% 
triple_backtrackFile = fopen('WT_NETseqTripleBacktrackParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==5)||(i==9)
            fprintf(triple_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*3);
        else
            fprintf(triple_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(triple_backtrackFile,'\n');
end
fclose(triple_backtrackFile);

%% 
quadruple_backtrackFile = fopen('WT_NETseqQuadrupleBacktrackParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==5)||(i==9)
            fprintf(quadruple_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*4);
        else
            fprintf(quadruple_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quadruple_backtrackFile,'\n');
end
fclose(quadruple_backtrackFile);

%% Make files of parameters to test
half_stallFile = fopen('WT_NETseqHalfStallParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==3)||(i==7)
            fprintf(half_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/2);
        else
            fprintf(half_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(half_stallFile,'\n');
end
fclose(half_stallFile);

%% 
third_stallFile = fopen('WT_NETseqThirdStallParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==3)||(i==7)
            fprintf(third_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/3);
        else
            fprintf(third_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(third_stallFile,'\n');
end
fclose(third_stallFile);

%% 
quarter_stallFile = fopen('WT_NETseqQuarterStallParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==3)||(i==7)
            fprintf(quarter_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/4);
        else
            fprintf(quarter_stallFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quarter_stallFile,'\n');
end
fclose(quarter_stallFile);

%% Make files of parameters to test
half_backtrackFile = fopen('WT_NETseqHalfBacktrackParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==5)||(i==9)
            fprintf(half_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/2);
        else
            fprintf(half_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(half_backtrackFile,'\n');
end
fclose(half_backtrackFile);

%% 
third_backtrackFile = fopen('WT_NETseqThirdBacktrackParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==5)||(i==9)
            fprintf(third_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)/3);
        else
            fprintf(third_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(third_backtrackFile,'\n');
end
fclose(third_backtrackFile);

%% 
quarter_backtrackFile = fopen('WT_NETseqQuarterBacktrackParams.txt','w');
for g=1:length(NETseq_geneNames)
    for i=1:size(NETseqGeneParamsBinned{model_idx},2)-3
        if (i==5)||(i==9)
            fprintf(quarter_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i)*0.25);
        else
            fprintf(quarter_backtrackFile,'%f\t',NETseqGeneParamsBinned{model_idx}(g,i));
        end
        
    end
    fprintf(quarter_backtrackFile,'\n');
end
fclose(quarter_backtrackFile);
