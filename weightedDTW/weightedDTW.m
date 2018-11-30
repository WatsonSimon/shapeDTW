% Dynamic Time Warping Algorithm
function [Dist,delta, D,op,rw,tw] = weightedDTW(t,r,pflag,metric,w)
%
% [Dist,D,k,w,rw,tw]=dtw(r,t,pflag,'absolute')
%
% Dist is unnormalized distance between t and r
% D is the accumulated distance matrix
% op: optimized path
% t: input path
% r: reference path
% pflag  plot flag: 1 (yes), 0(no)
% w: window parameter, if s(i) is matched with t(j) then |i-j|<=w
% rw is the warped r vector
% tw is the warped t vector
% metric: method to calculate distance
% 0: absolute
% 1: square

%% Normalize 2 inputs
% r = zscore(r);
% t = zscore(t);

%% check the size of input1
[row,M] = size(r); 
if (row > M) 
    M = row; 
    r = r'; 
end;

%% check the size of input2
[row,N] = size(t); 
if (row > N) 
    N = row; 
    t = t'; 
end;

%% Calculate pair-wise distance between input1 and input2

% Calculate weighted vector
length = max(size(t,2), size(r,2));
for i = 1:length+1
    wei(i) = 1/(1+exp(-0.1*((i-1)-length/2)));
end

d = 0;
for i=1:row
    tt = t(i,:);
    rr = r(i,:);
    tt = (tt-mean(tt))/std(tt);
    rr = (rr-mean(rr))/std(rr);
    
%     if metric == 0
% %         d = abs(repmat(r',1,N)-repmat(t,M,1));   
%         d = d + abs(repmat(rr',1,N)-repmat(tt,M,1)); 
%     elseif metric == 1
% %         d = (repmat(r',1,N)-repmat(t,M,1)).^2;
%         d = d + (repmat(rr',1,N)-repmat(tt,M,1)).^2;
%     end

    for i = 1:size(rr,2)
        for j = 1:size(tt,2)
            ii = abs(i-j) + 1;
            d(i,j) = wei(ii)*(rr(i) - tt(j))^2;
        end
    end
end

%% Initialze a new distance matrix
D = zeros(size(d))+Inf;
D(1,1) = 0;

%% Update table D(i;j): Accumulated cost
% 2 special upper border cases 
% Other cases
for m = 2:M
    for n = max(m-w,2):min(m+w,N)
        D(m,n) = d(m,n) + min([D(m-1,n),D(m-1,n-1),D(m,n-1)]);
    end
end
 
Dist = D(M,N);
n = N; m = M; k = 1;
op = [M N];

%% Dynamic programming
while ((n+m)~=2) % Until (n+m) = 2
    if n == 1
        m = m-1;
    elseif m == 1
        n = n-1;
    else 
      [values,index] = min([D(m-1,n),D(m,n-1),D(m-1,n-1)]);
      
      switch index
        case 1
            m = m-1;
        case 2
            n = n-1;
        case 3
            m = m-1;
            n = n-1;
      end
  end
    k = k+1;
    op = [m n; op]; % Add index of matched pairs
end
 
%% warped path
rw = r(op(:,1));
tw = t(op(:,2));

%% Calculate summative score
delta = 0;
alpha = 1.019;
beta = 1.171;
gamma = 1.553;

for i = 1:size(rw,2)
    d1 = abs(r(op(i,1)) - t(op(i,2)));
    k1 = 1;
    delta1 = (alpha^((gamma^-d1)*k1))*(beta^d1);
    delta = delta + delta1;
end
 
%% Choose to draw the result
if pflag
    
    % --- Accumulated distance matrix and optimal path
    figure('Name','DTW - Accumulated distance matrix and optimal path', 'NumberTitle','off');
    
    main1 = subplot('position',[0.19 0.19 0.67 0.79]);
    image(D);
    cmap = contrast(D);
    colormap(cmap); % 'copper' 'bone', 'gray' imagesc(D);
    hold on;
    x = op(:,1); y = op(:,2);
    ind=find(x==1); x(ind)=1+0.2;
    ind=find(x==M); x(ind)=M-0.2;
    ind=find(y==1); y(ind)=1+0.2;
    ind=find(y==N); y(ind)=N-0.2;
    plot(y,x,'-w', 'LineWidth',1);
    hold off;
    axis([1 N 1 M]);
    set(main1, 'FontSize',7, 'XTickLabel','', 'YTickLabel','');
 
    colorb1 = subplot('position',[0.88 0.19 0.05 0.79]);
    nticks=8;
    ticks=floor(1:(size(cmap,1)-1)/(nticks-1):size(cmap,1));
    mx=max(max(D));
    mn=min(min(D));
    ticklabels = floor(mn:(mx-mn)/(nticks-1):mx);
    colorbar(colorb1);
    set(colorb1, 'FontSize',7, 'YTick',ticks, 'YTickLabel',ticklabels);
    set(get(colorb1,'YLabel'), 'String','Distance', 'Rotation',-90, 'FontSize',7, 'VerticalAlignment','bottom');
    % Draw signal t
    left1 = subplot('position',[0.07 0.19 0.10 0.79]);
    plot(r,M:-1:1,'-b');
    set(left1, 'YTick',mod(M,10):10:M, 'YTickLabel',10*rem(M,10):-10:0)
    axis([min(r) 1.1*max(r) 1 M]);
    set(left1, 'FontSize',7);
    set(get(left1,'YLabel'), 'String','Samples', 'FontSize',7, 'Rotation',-90, 'VerticalAlignment','cap');
    set(get(left1,'XLabel'), 'String','Amp', 'FontSize',6, 'VerticalAlignment','cap');
    % Draw signal r
    bottom1 = subplot('position',[0.19 0.07 0.67 0.10]);
    plot(t,'-r');
    axis([1 N min(t) 1.1*max(t)]);
    set(bottom1, 'FontSize',7, 'YAxisLocation','right');
    set(get(bottom1,'XLabel'), 'String','Samples', 'FontSize',7, 'VerticalAlignment','middle');
    set(get(bottom1,'YLabel'), 'String','Amp', 'Rotation',-90, 'FontSize',6, 'VerticalAlignment','bottom');
    
    %% --- Warped signals
    figure('Name','DTW - warped signals', 'NumberTitle','off');
    
    subplot(1,2,1);
    set(gca, 'FontSize',7);
    hold on;   
    plot(r,'-bx'); plot(t,':r.');
    hold off;
    axis([1 max(M,N) min(min(r),min(t)) 1.1*max(max(r),max(t))]);
    grid;legend('TS1','TS2'); title('Original'); xlabel('Time'); ylabel('Value');
    
    subplot(1,2,2);
    set(gca, 'FontSize',7);
    hold on;
    plot(rw,'-bx'); plot(tw,':r.');
    hold off;
    axis([1 k min(min([rw; tw])) 1.1*max(max([rw; tw]))]);
    grid;legend('TS1','TS2'); title('Optimum'); xlabel('Time'); ylabel('Value'); 
end