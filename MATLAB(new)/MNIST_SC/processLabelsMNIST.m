function X = processLabelsMNIST(filename)


dataFolder = fullfile(tempdir,'mnist');
gunzip(filename,dataFolder)

[~,name,~] = fileparts(filename);

[fileID,errmsg] = fopen(fullfile(dataFolder,name),'r','b');
if fileID < 0
    error(errmsg);
end

magicNum = fread(fileID,1,'int32',0,'b');
if magicNum == 2049
    fprintf('\nRead MNIST label data...\n')
end

numImages = fread(fileID,1,'int32',0,'b');
fprintf('Number of labels in the dataset: %6d ...\n',numImages);

X = fread(fileID,inf,'unsigned char');

X = reshape(X,[1,size(X,1)]);

%X = reshape(X,numCols,numRows,numImages);
%X = permute(X,[2 1 3]);
%X = X./255;
%X = reshape(X, [28,28,1,size(X,3)]);

fclose(fileID);
end