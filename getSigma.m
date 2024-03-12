function force = getSigma(force,px2m)
    x = cat(1,force.x);y=cat(1,force.y);z=cat(1,force.z);r=cat(1,force.r);rr=r/px2m;
    for i=1:numel(force)
        i;
%         %%this was for testing the rotation of beta.
%             r=r/px2m;
%             if z(i)>0
%                 figure(1);clf(1);hold on;
%                 n = force(i).n;
%                 rectangle('Position',[x(i)-rr(i),y(i)-rr(i),2*rr(i),2*rr(i)],'EdgeColor','r','Curvature',[1,1]);
%                 col = colormap(hsv(z(i)));
%                 b=force(i).beta+pi;
%                 for j=1:z(i)
%                     jj=n(j);
%                     if jj>0
%                     rectangle('Position',[x(jj)-rr(jj),y(jj)-rr(jj),2*rr(jj),2*rr(jj)],'EdgeColor',col(j,:),'Curvature',[1,1]);            
%                     plot([x(i) x(i)+rr(jj)*cos(b(j))],[y(i) y(i)+rr(jj)*sin(b(j))],'-','color',col(j,:));
%                     end
%                 end
%                 pause
%             end                  
        if (z(i)>0) & ~isempty(force(i).forces)
            n = force(i).neighbours;
            rij = zeros(z(i),2);
            for j=1:z(i)
                if (numel(force)>n(j)&&n(j)>0)
                   
                    rij(j,1)  = (x(n(j))-x(i))*px2m/2;
                    rij(j,2)  = (y(n(j))-y(i))*px2m/2;
                else
                    rij(j,1)  = force(i).r*cos(force(i).betas(j)+pi);
                    rij(j,2)  = force(i).r*sin(force(i).betas(j)+pi);                    
                end        
            end  
            f = force(i).forces; alpha = force(i).alphas; beta = force(i).betas+pi;
            fx = f.*cos(alpha+beta);
            fy = f.*sin(alpha+beta);
            fxy = [fx',fy'];
            %--- stress tensor            
            M = zeros(2,2,z(i));
            for j=1:z(i)
                
                M(:,:,j) =  rij(j,:)'*fxy(j,:);
            end            
            sigma = sum(M,3);
            force(i).sigma = sigma;                     
        else
            sigma = nan(2,2);
            force(i).sigma = sigma;
        end
    end
end

