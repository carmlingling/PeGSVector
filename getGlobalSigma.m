function [Sigma] = getGlobalSigma(force,px2m)
Sigma = zeros(2,2);
%     for i = 1:numel(force)
%         if ~isnan(force(i).sigma)
%             Sigma = Sigma + force(i).sigma;
%         end
%     end

% particle position list
x = cat(1,force.x); y=cat(1,force.y); z=cat(1,force.z);

% look at all particles
for i=1:numel(force)
if (z(i)>0) & ~isempty(force(i).forces)
        n = force(i).neighbours;
        % look at all neighbors with larger indices
        for j=1:numel(n)
            rij = zeros(1,2);
            if n(j) > i
                % center-to-center vector
                rij(1,1)  = (x(n(j))-x(i))*px2m/2;
                rij(1,2)  = (y(n(j))-y(i))*px2m/2;
                % force vector
                f = force(i).forces(j); alpha = force(i).alphas(j); beta = force(i).betas(j)+pi;
                fx = f.*cos(alpha+beta);
                fy = f.*sin(alpha+beta);
                fxy = [fx',fy'];
                
                Sigma = Sigma + real(rij'*fxy);
            end
        end
    end
end
end
