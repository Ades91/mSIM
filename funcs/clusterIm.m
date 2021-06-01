% takes a binary image as input
% return a list of clustered pixels
function [out,clusteredImage] = clusterIm(input,clusterSizeThreshold, display)

if nargin < 3; display = 0; end
if nargin < 2; clusterSizeThreshold = numel(input); end

clusterList = [];
tempClusterList = [];
visited = zeros(size(input)); % keep track of unvisited pixels
neighList = [];
cPixel = [];
for yy = 1:size(input,1)
    for xx = 1:size(input,2)
        if (input(yy,xx) ==1 && visited(yy,xx) == 0) % if we find and unvisited pixel
            clusterList{end+1} = [xx,yy]; visited(yy,xx) = 1; % create a new cluster and mark as visited
            % check is we can find unvisited neighbourg 
            neighList = [neighList; getNeighbourg(xx,yy,input,visited)];
            while ~isempty(neighList)
                cPixel = neighList(1,:); neighList(1,:) = []; % pick a pixel from the list
                if (visited(cPixel(2),cPixel(1)) == 0)
                    clusterList{end}(end+1,:) = cPixel; visited(cPixel(2),cPixel(1)) = 1; %attach to cluster and mark as visited
                    neighList = [neighList; getNeighbourg(cPixel(1),cPixel(2),input,visited)]; % update neighList using cPixel pos
                end
            end
        end
    end
end

% get cluster sizes
for k = 1:length(clusterList)
    cSize(k) = length(clusterList{k});
end
[~,map] = sort(cSize,'descend');

if clusterSizeThreshold < 1
    % remove a fraction of the smallest cluster based on clusterSizeThreshold
    map = map(1:round(clusterSizeThreshold*length(clusterList))); %keep indices of top x% biggest cluster
    out = cell(length(map),1);
    for k = 1:length(map)
        out{k} = clusterList{map(k)}; 
    end
else
    nCluster = sum(cSize >= clusterSizeThreshold);
    out = cell(nCluster,1);
    for k = 1:nCluster
        out{k} = clusterList{map(k)}; % keep the nCluster larger than clusterSizeThreshold
    end
end
    
clusteredImage = zeros(size(input));
for k = 1:length(out)
    for j = 1:size(out{k},1)
        clusteredImage(out{k}(j,2),out{k}(j,1)) = k;
    end
end


if display
 

   figure(display); imagesc(clusteredImage)
end

end

function nei = getNeighbourg(x,y,input,visited)
    nei = [];
    
    list = [x,y-1; x,y+1; x-1,y; x+1,y];                    % 4-connected
%     list = [list; x-1,y-1; x-1,y+1; x+1,y-1; x+1,y+1];    % 8-connected
        
    for k = 1:size(list,1)
        if list(k,2) > 0 && list(k,2) <=  size(input,1) && list(k,1) > 0 && list(k,1) <= size(input,2)
            try
            if (visited(list(k,2),list(k,1))==0 && input(list(k,2),list(k,1)) == 1)
                nei(end+1,:) = list(k,:);
            end
            catch
               disp('oups') 
            end
        end
    end
%     
%     % bottom pixel
%     if y < length(input)
%         if (visited(y+1,x)==0 && input(y+1,x) == 1)
%             nei(end+1,:) = [x,y+1];
%         end
%     end
%     % left pixel
%     if x > 1
%         if (visited(y,x-1)==0 && input(y,x-1) == 1)
%             nei(end+1,:) = [x-1,y];
%         end
%     end
%     % right pixel
%     if x < length(input)
%         if (visited(y,x+1)==0 && input(y,x+1) == 1)
%             nei(end+1,:) = [x+1,y];
%         end
%     end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% DEPRECATED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deprecated way to determine if a pixel belongs to an existing cluster
%             c = isConnected(tempClusterList,xx,yy);
%            if c > 0
%                tempClusterList{c}(end+1,:) = [xx,yy];
%            else
%                tempClusterList{end+1} = [xx,yy];
%            end
function c = isConnected(cluster,x,y)
c = 0;
for k = 1:length(cluster)
%     disp(['Cluster ', num2str(k),' ************'])
    clus = cluster{k};
   for j = 1:size(clus,1)
       if  (((clus(j,1)-x)^2 + (clus(j,2)-y)^2)) <= 3
          c = k;
          break;
       end
   end
end

end
%% scnd phase : consolidate the clusters
% 
% count = 0;
%     currentCluster = tempClusterList{1};
%     tempClusterList{1} = []; tempClusterList(~cellfun('isempty',tempClusterList));
%     
% while (~isempty(tempClusterList) && count < 1000)
% 
%     for k = 1:length(tempClusterList)
%         for currentInd = 1:size(currentCluster,1)
%            for tempInd = 1:size(tempClusterList{k},1)
%                if ((currentCluster(currentInd,1)-tempClusterList{k}(tempInd,1))^2 + ...
%                        (currentCluster(currentInd,1)-tempClusterList{k}(tempInd,1)^2)) <= 3
%                    
%                    
%                    
%                end
%            end
%         end
%     end
%     count = count+1;
% end