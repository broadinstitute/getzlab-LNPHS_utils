function C = interval_list_intersect(A,B)
% C = interval_list_intersect(A,B)
%
% give two interval lists A and B with format [start1 end1; start2 end2; ... ; startN endN],
% computes their intersection C with the same format
%
% Mike Lawrence 2010-10-05

if nargin~=2, error('takes two args: A and B'); end
if size(A,2)~=2 || size(B,2)~=2, error('A and B should have two columns each'); end
if size(A,1)==0 || size(B,1)==0, C = zeros(0,2); return; end

x = [A(:);B(:)];
mn = min(x);
mx = max(x);
ln = mx-mn+1; 
adj = mn-1;

% current implementation is a memory hog but fast

d = zeros(ln,1,'uint8');
for i=1:size(A,1)
  if A(i,1)>A(:,2), error('start>end'); end
  d(A(i,1)-adj:A(i,2)-adj) = bitor(1,d(A(i,1)-adj:A(i,2)-adj));
end
for i=1:size(B,1)
  if B(i,1)>B(:,2), error('start>end'); end
  d(B(i,1)-adj:B(i,2)-adj) = bitor(2,d(B(i,1)-adj:B(i,2)-adj));
end

d = 1*(d==3);
b = diff([0;d;0]);
C = [find(b==1) find(b==-1)-1] + adj;


