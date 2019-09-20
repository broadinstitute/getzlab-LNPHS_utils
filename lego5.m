function lego5(z,P)
% input should be a (8x4)x(12x4) matrix of heights for the "LEGO5" plot
% tries to convert other formats
%
% blank spaces for alignment with lego.m
%

if ~exist('P','var'), P=[]; end

if iscellstr(P)   % can specify catnames as second argument for factor2lego
  tmp=P;
  P=[];
  P.catnames = tmp;
end

if iscell(z), lego5s(z,P); return; end
if isstruct(z) && isfield(z,'mut'), z=z.mut; end
if ischar(z) || (isstruct(z) && (isfield(z,'context1025')||(isfield(z,'chr')&(isfield(z,'st')|isfield(z,'start')|isfield(z,'pos'))))), lego5maf(z,P); return; end
%
% blank spaces for alignment with lego.m
%
if ~(length(size(z))==2 & all(size(z)==[8*4 12*4])), error('unknown format'); end
 
if ischar(P) && ismember(P,{'counts','rates','exome','genome'})
  if ~strcmp(P,'counts'), error('only "counts" supported'); end
  P=[];
end

P = impose_default_value(P,'log',false);
P = impose_default_value(P,'imagesc',false);  % if true, displays as imagesc instead of bar3_with_colors
P = impose_default_value(P,'zmax',[]);

if P.log, z=log(z); end

if P.imagesc
  if ~isempty(P.zmax), z = min(z,P.zmax); end
  imagesc(z);
  pp = {'linewidth',0.5,'color',[0 0 0]};
  for x=4.5:4:47.5, line([x x],ylim,pp{:}); end
  for y=4.5:4:31.5, line(xlim,[y y],pp{:}); end
  pp = {'linewidth',3,'color',[0 0 0]};
  line([16.5 16.5],ylim,pp{:});line([32.5 32.5],ylim,pp{:});line(xlim,[16.5 16.5],pp{:});
else
  c = get_LEGO5_colors();
%  z3 = repmat(z,[1 1 3]); c(isnan(z3))=1; % white   %%% for some reason this is causing bizarre color scrambling in context of legos.m
  bar3_with_colors(z,1,c);
  if ~isempty(P.zmax), zlim([0 P.zmax]); end
  % NOTE: sometimes in Mac Xwin "Xquartz", the plot disappears if zmax drops below ~0.4e-6
  zl = zlim;
  if zl(2)<0.4e-6
    fprintf('Warning: plot sometimes disappears when zmax<0.4e-6.  If this happens, try typing zlim([0 0.4e-6]).\n');
  end
end
set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);



