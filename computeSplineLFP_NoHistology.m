function [splineLFP, xColumns,soma_depth] = computeSplineLFP_NoHistology(mnLFP, mnLFP2, equalChans, info, winAroundSpike) % difference between mnLFP and mnLFP2?

test = size(info.connected,1);
    mnLFP = mnLFP(info.connected,:); % first getting rid of reference channels (flat signal) or broken channels
    mnLFP2 = mnLFP2(info.connected,:);
    
    % standardising
    
%     mnLFP = mnLFP ./ std_vector_TESTING';
%     mnLFP2 = mnLFP2 ./ std_vector_TESTING';

    
    % to test if they were succesfully removed you can test:
%     figure
%     plot(std(mnLFP'))
      
    spikeTime = (size(winAroundSpike,2) +1)/2; % this is similar to soma time
    [~, chanMin] = min(mnLFP(:,spikeTime)); % this is soma channel
    soma_depth = chanMin;
    xExtent = [max(info.chanCoords(:,1)), min(info.chanCoords(:,1))];
    shift = [-16,16]; 

    % next section determines what column the smallest amplitude is:
            
            % channel map:
            
            %           1 2 3 4         (columns)
            %
            %          |0   0  | 
            %          |  0   0|
            %          |0   0  |
            %          |  0   0|            
            %          |0   0  |
            %          |  0   X|        (X - example channel with  highest amplitude)
            %          |0   0  |
            %          |  0   0|           
            %          |0   0  |
            %          |  0   0|        % also important: reference chans!
            %          \       /
            %           \     /
            %            \   /
            %             \ /
            %              ▼
            
            % in the case example above (X), the first if clause is run
            % info.chanCoords(chanMin,1) gives the coordinate 
        
        
    if info.chanCoords(chanMin,1) == xExtent(1) % if highest amplitude is on the farthest right column (column 4, for example where the X is)
        ys = info.chanCoords(info.chanCoords(:,1) == info.chanCoords(chanMin,1),2); % (column 4 y coords) -
        ys2 =  info.chanCoords(info.chanCoords(:,1) == (info.chanCoords(chanMin,1) + shift(1)),2); % (column 3 y)
        id = info.chanCoords(chanMin,1); % takes the column with max amplitude
        id2 = info.chanCoords(chanMin,1) + shift(1); % and the one next to it (in this case +16)
        
        
        
    elseif info.chanCoords(chanMin,1) == xExtent(2) % if its on the other extreme (column 1) 
        ys = info.chanCoords(info.chanCoords(:,1) == info.chanCoords(chanMin,1),2); % (column 1 y coords)
        ys2 =  info.chanCoords(info.chanCoords(:,1) == (info.chanCoords(chanMin,1) + shift(2)),2); % (column 2 y coords)
        id = info.chanCoords(chanMin,1); % takes the column with max amplitude - this is only a scalar - x coordinates have same
        id2 = info.chanCoords(chanMin,1) + shift(2); % and the one next to it (in this case -16)
        
        
        
	else  % get min value at columns nearest to where the soma is (chanMin)
       [valMin2, ~] = min(mnLFP(info.chanCoords(chanMin,1)-16,spikeTime)); % checks the max amplitude on two neighbouring columns
       [valMin3, ~] = min(mnLFP(info.chanCoords(chanMin,1)+16,spikeTime));
         
       [~, idx] = min([valMin2,valMin3]); %smallest of the two neighbouring channels
       ys = info.chanCoords(info.chanCoords(:,1) == info.chanCoords(chanMin,1),2); % the one with max
       ys2 =  info.chanCoords(info.chanCoords(:,1) == (info.chanCoords(chanMin,1)+ shift(idx)),2);% determines the shift based on idx (which adjecent channel is  chosen - second highest amplitude
       id = info.chanCoords(chanMin,1);
       id2 = info.chanCoords(chanMin,1) + shift(idx); % same as for ys2
    end

    ysN = nearestpoint(min(ys):20:max(ys),equalChans);
     newysN = nearestpoint(min(ys2):20:max(ys2),equalChans);     
    
    
    xColumns = [id id2]; % just to be on the safe side that it is actually taking the same columbns for different spike numbers
    
  
    % make variable here that will include all ys - even if column is
    % switched, smoething like min(ys):20:max(ys)
    
    for j = 1:size(mnLFP2,2) % this now defines voltage values in each column as the ones from the two adjecent columns
%       try

	sM = NaN(2,size(min(equalChans):20:max(equalChans),2));
      sM(1,ysN) = csapi(ys, mnLFP2(info.chanCoords(:,1) == id,j), min(ys):20:max(ys)); % qubic spline inrepolation along distance for first column
      sM(2,newysN) = csapi(ys2, mnLFP2(info.chanCoords(:,1) == id2,j), min(ys2):20:max(ys2)); % along second column

        %NOTE : this section can be used to increase spatial respolution by
        % interpolating
        
%               sM = NaN(2,size(min(ys):2:max(ys),2));
%         sM(1,:) = csapi(ys, mnLFP2(info.chanCoords(:,1) == id,j), min(ys):2:max(ys)); % qubic spline inrepolation along distance for first column
%         sM(2,:) = csapi(ys2, mnLFP2(info.chanCoords(:,1) == id2,j), min(ys2):2:max(ys2)); % along second column
%     
      splineLFP(:,j) = nanmean(sM,1).';
% 	  catch
% 		  continue
% 	  end

    end
    
    
 end