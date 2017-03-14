function D = Produce_data(c, X1)
% produce data table randomly from training image
j = 0;
for i=1:size(X1,2)
    data1((j+1):(i*size(X1,1)),1) = 1:size(X1,1); % Y co-ordinate of data points
    data1((j+1):(i*size(X1,1)),2) = i; % X co-ordinate of data points
    data1((j+1):(i*size(X1,1)),3) = X1(:,i); % Conditioning values
    j = j + size(X1,1);
end
data = datasample(data1,c,'Replace',false);

% data1 = datasample(X,1,2,'Replace',false);
% data2 = datasample(X,1,1,'Replace',false);
% data = data1;
% data((size(data1,1)+1):(size(data1,1)+size(data2,2))) = data2';
% data = datasample(data,c,'Replace',false);
% % xlocation = randi(size(X,1),c,1);
% % ylocation = randi(size(X,2),c,1);

% xlocation = randi([2,p-1],c,1);
% ylocation = randi([2,q-1],c,1);
% D(:,1) = xlocation;
% D(:,2) = ylocation;
% D(:,3) = data;

D = data;