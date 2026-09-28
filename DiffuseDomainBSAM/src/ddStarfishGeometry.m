function geom = ddStarfishGeometry(n)
%DDSTARFISHGEOMETRY The five-armed starfish in DiffuseDomainHelmholtz.
% r(theta) = 0.9*(1.2 + 0.7*sin(5*theta)). Default: 3000 boundary points,
% including the repeated endpoint, matching the reference generateBoundary.
% Distance is the exact point-to-segment distance to this polygon (up to
% rounding), positive inside. A bounding-box tree prunes distant segments.
% This needs neither a MEX build nor a Statistics/Image Processing toolbox.
if nargin < 1, n = 3000; end
validateattributes(n, {'numeric'}, {'scalar','integer','>=',32});
theta = linspace(0, 2*pi, n).';
r = 0.9*(1.2 + 0.7*sin(5*theta));
v = [r.*cos(theta), r.*sin(theta)]; v(end,:) = v(1,:);
a = v(1:end-1,:); b = v(2:end,:); e = b-a;
len2 = sum(e.^2,2); ns = n-1;
nodes = struct('box',{},'left',{},'right',{},'first',{},'last',{});
build(1, ns);
geom.name = 'Helmholtz starfish: r = 0.9(1.2 + 0.7 sin(5 theta))';
geom.vertices = v;
geom.signedDistance = @distance;
geom.psi = @(X,Y) -distance(X,Y);
geom.inside = @(X,Y) distance(X,Y) >= 0;

    function id = build(first,last)
        id = numel(nodes)+1;
        vv = [a(first:last,:); b(first:last,:)];
        nodes(id) = struct('box',[min(vv,[],1),max(vv,[],1)], ...
            'left',0,'right',0,'first',first,'last',last);
        if last-first > 7
            mid = floor((first+last)/2);
            left = build(first,mid); right = build(mid+1,last);
            nodes(id).left = left; nodes(id).right = right;
        end
    end

    function d = distance(X,Y)
        assert(isequal(size(X),size(Y)), 'X and Y must have equal size.');
        shape = size(X); x = X(:); y = Y(:);
        angle = mod(atan2(y,x),2*pi);
        seg = min(floor(angle*ns/(2*pi))+1, ns);
        % A boundary vertex gives a valid initial distance upper bound.
        best = (x-a(seg,1)).^2 + (y-a(seg,2)).^2;
        stackNode = {1}; stackIdx = {(1:numel(x)).'};
        while ~isempty(stackNode)
            id = stackNode{end}; stackNode(end) = [];
            idx = stackIdx{end}; stackIdx(end) = [];
            node = nodes(id); bb = node.box;
            dx = max(max(bb(1)-x(idx),0),x(idx)-bb(3));
            dy = max(max(bb(2)-y(idx),0),y(idx)-bb(4));
            idx = idx(dx.^2+dy.^2 <= best(idx)+64*eps);
            if isempty(idx), continue; end
            if node.left ~= 0
                stackNode(end+1:end+2) = {node.left,node.right};
                stackIdx(end+1:end+2) = {idx,idx};
            else
                xx=x(idx); yy=y(idx); local=best(idx);
                for j=node.first:node.last
                    t=((xx-a(j,1))*e(j,1)+(yy-a(j,2))*e(j,2))/len2(j);
                    t=max(0,min(1,t));
                    dd=(xx-a(j,1)-t*e(j,1)).^2+(yy-a(j,2)-t*e(j,2)).^2;
                    local=min(local,dd);
                end
                best(idx)=local;
            end
        end
        % This polygon is radially ordered. Intersect the ray at theta with
        % its enclosing segment to classify the polygon interior exactly.
        ux=cos(angle); uy=sin(angle);
        radius=(a(seg,1).*e(seg,2)-a(seg,2).*e(seg,1)) ./ ...
            (ux.*e(seg,2)-uy.*e(seg,1));
        s=ones(size(x)); s(hypot(x,y)>radius)=-1;
        d=reshape(s.*sqrt(best),shape);
    end
end
