%particle emission:

%2 populations: reflected and sputtered particles
%a) sputtered are emitted with E=Thompson distribution &
%   angle=under/over-cosine (to test sensitivity)
%b) reflected are emitted with E=energy reflection coeff. & angle=near
%   specular reflection

%particles travel in straight lines, so we can use functions built for
%launching initial particles here.

fprintf(1, '\n');
disp('---- calling particle emission -----')

disp('re-deposition:   0:  YES, particle re-deposited in final position')
disp('                1-6: NO; Final position = projection of trajectory in initial z0: ')
disp('                1:  (x0,y0)=(xs,yx)     2: err in z-direction  3: err in x- or y-direction;')
disp('                5: invalid exitflag    6:  fzero(Exitflag)=6     7: invalid theta  ')



zk(1:npoints,1:npoints)=0.0;

%the function that finds the intersection with the surface is similar to
%that for the initial particles, but here we need also need to look for
%particles travelling in +z, which might hit the surface, or not

nSppart(1:npoints,1:npoints)=0.0; %number of particles emitted by sputtering
SpEout(1:npoints,1:npoints,:)=0.0; % energy of particles emitted through sputtering
SpThout(1:npoints,1:npoints,:)=0.0; %theta angle of particles emitted by sputtering
SpThout_loc(1:npoints,1:npoints,:)=0.0;
aveSpThout(1:npoints,1:npoints)=0.0;
aveSpThout_loc(1:npoints,1:npoints)=0.0;
SpPhiout_loc(1:npoints,1:npoints,:)=0.0;
SpPhiout(1:npoints,1:npoints,:)=0.0; %phi angle of particles emitted by sputtering
aveSpPhiout(1:npoints,1:npoints)=0.0;
aveSpPhiout_loc(1:npoints,1:npoints)=0.0;
sput_redepos_p(1:4)=0.0;
sput_redepos(1:npoints,1:npoints,:,1:4)=0.0;


%idem for particles emitted through reflection
nRpart(1:npoints,1:npoints)=0.0;
REout(1:npoints,1:npoints)=0.0;
RThout(1:npoints,1:npoints,:)=0.0; %theta angle of particles emitted by sputtering
RThout_loc(1:npoints,1:npoints,:)=0.0;
aveRThout(1:npoints,1:npoints)=0.0;
aveRThout_loc(1:npoints,1:npoints)=0.0;
RPhiout_loc(1:npoints,1:npoints)=0.0;
RPhiout(1:npoints,1:npoints)=0.0; %phi angle of particles emitted by reflection
refl_redepos_p(1:4)=0.0;
refl_redepos(1:npoints,1:npoints,:,1:4)=0.0;
NredepSp(1:npoints,1:npoints)=0.0; %where sputtered particles are redeposited
NredepR(1:npoints,1:npoints)=0.0;  %where sputtered particles are redeposited


for i=1:npoints
    for j=1:npoints
        
        %number of particles emitted from (i,j) cell that are redeposited
        NSp=0;
        NR=0;
        
        %%%%%%%%%%% 1-Cell characteristics %%%%%%%%%%%
        xi=xg(i)+0.5*dx;
        yj=yg(j)+0.5*dy;
        zk(i,j)=zs(xi,yj,Ax,fx,Ay,fy);
        
        
        %%%%%%%%%%%%%%%% 2-Sputtering %%%%%%%%%%%%%%%%
        nSppart(i,j)=Y_loc(i,j)*Ncounts(i,j);
        rnSp=rand;
        if (( nSppart(i,j) - floor(Y_loc(i,j)*Ncounts(i,j)) ) <= rnSp )
            nEmittSp=floor(nSppart(i,j));            
        elseif (( nSppart(i,j) - floor(Y_loc(i,j)*Ncounts(i,j)) ) > rnSp )
            nEmittSp=floor(nSppart(i,j))+1;
        else
            disp('ERR nEmittSp')
            disp(num2str(nSppart(i,j)))
        end
            
        nSppart(i,j)=nEmittSp;

        
        % energy and angle of particles emitted through sputtering, wrt
        % local normal to surf
        % all reflected particles are emitted with same energy, REout(i,j)
        if(nEmittSp> 0.0)
            
            %how many particles emitted from this cell were re-deposited
            NSp=0.0;
            NR=0.0;
            
            %Stest_nemitt=['for cell in x=', num2str(xi), ' y=',
            %num2str(yj) , ' ; nEmittR=', num2str(nEmittSp), ' (RN_loc=',num2str(Y_loc(i,j)) , ' , Ncounts=',num2str(Ncounts(i,j)), ')' ];
            %disp(Stest_nemitt)
            
            %initialize
            SpThout(i,j,1:nEmittSp)=0.0; %theta angle of particles emitted by sputtering
            SpThout_loc(i,j,1:nEmittSp)=0.0;
            SpPhiout_loc(i,j,1:nEmittSp)=0.0;
            SpPhiout(i,j,1:nEmittSp)=0.0; %phi angle of particles emitted by sputtering
            sput_redepos(i,j,1:nEmittSp,1:4)=0.0;
            aveSpThout_loc_tmp=0.0;
            aveSpThout_tmp=0.0;
            aveSpPhiout_loc_tmp=0.0;
            aveSpPhiout_tmp=0.0;
            
            SpEout(i,j,1:nEmittSp) = thompson_distr(E0, Esb, Mtg, Mpr, nEmittSp, npoints); %needed?
            SpThout_loc(i,j,1:nEmittSp) = cosn_distrib(r1, n1, r2, n2, nEmittSp, npoints);
            %disp(180/pi*SpThout_loc(i,j,1:nEmittSp))
            
            %For each sputtered particle
            for p=1:nEmittSp
                %define characteristics, wrt global coordinate system
                SpThout(i,j,p)=sg_theta(i,j)+SpThout_loc(i,j,p);                
                
                
                if (sign(sg_theta(i,j))==sign(SpThout_loc(i,j,p))) %adding angles -> random, within (surface phi +/- pi/2)
                    SpPhiout_loc(i,j,p)=(rand-0.5)*pi; %
                else %resting angles for SpThout -> random, within (surface phi +/- pi)
                    SpPhiout_loc(i,j,p)=(rand-0.5)*2*pi;
                end
                SpPhiout(i,j,p)=sg_phi(i,j)+SpPhiout_loc(i,j,p);
                
                
                %to take averages
                aveSpThout_loc_tmp=aveSpThout_loc_tmp+SpThout_loc(i,j,p);
                aveSpThout_tmp=aveSpThout_tmp+SpThout(i,j,p);
                aveSpPhiout_loc_tmp=aveSpPhiout_loc_tmp+SpPhiout_loc(i,j,p);
                aveSpPhiout_tmp=aveSpPhiout_tmp+SpPhiout(i,j,p);
                
                %intersection of sputtered particle with surface
                sput_redepos_p=emitted_part_surf_intersec(xi,yj,zk(i,j),SpPhiout(i,j),SpThout(i,j,p),Ax,fx,Ay,fy,z0,dx, dy);
                sput_redepos(i,j,p,1:4)=sput_redepos_p(1:4);   %to store the data
                xr=sput_redepos_p(1);   %for sorting redep positions
                yr=sput_redepos_p(2);
                zr=sput_redepos_p(3);
                kr=sput_redepos_p(4);
                
                
                %{
                S7=['     particle sputtered from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with theta =', num2str(180/pi*SpThout(i,j,p)), ' and phi =', num2str(180/pi*SpPhiout(i,j,p))];
                S8=['     re-deposited ', num2str(sput_redepos(i,j,p,4)) , ' ; Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                disp(S7)
                disp(S8)
                fprintf(1, '\n');
                %}
                
                if (kr==0)
                    
                    ir=floor((xr-surfxmin)/di)+1;
                    if (ir<=0)
                        %{
                        errir=['WARNING, redeposited out of surface area; xr=', num2str(xr), ' ir=', num2str(ir) ' ; setting ir=1' ];
                        disp(errir)
                        S7=['     particle sputtered from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*SpThout_loc(i,j,p)) ' theta =', num2str(180/pi*SpThout(i,j,p)), ' local phi = ' num2str(180/pi*SpPhiout_loc(i,j,p)) ' and phi =', num2str(180/pi*SpPhiout(i,j))];
                        S8=['     re-deposited ', num2str(sput_redepos(i,j,p,4)) , ' ; Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        ir=1;
                    elseif (ir>npoints)
                        %{
                        errir=['WARNING, redeposited out of surface area; xr=', num2str(xr), ' ir=', num2str(ir) ' ; setting ir=npoints' ];
                        disp(errir)
                        S7=['     particle sputtered from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*SpThout_loc(i,j,p)) ' theta =', num2str(180/pi*SpThout(i,j,p)), ' local phi = ' num2str(180/pi*SpPhiout_loc(i,j,p)) ' and phi =', num2str(180/pi*SpPhiout(i,j))];
                        S8=['     re-deposited ', num2str(sput_redepos(i,j,p,4)) , ' ; Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        ir=npoints;
                    end
                    
                    jr=floor((yr-surfymin)/dj)+1;
                    if (jr<=0)
                        %{
                        errjr=['WARNING, redeposited out of surface area; yr=', num2str(yr), ' jr=', num2str(jr) ' ; setting jr=1' ];
                        disp(errjr)
                        S7=['     particle sputtered from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*SpThout_loc(i,j,p)) ' theta =', num2str(180/pi*SpThout(i,j,p)), ' local phi = ' num2str(180/pi*SpPhiout_loc(i,j,p)) ' and phi =', num2str(180/pi*SpPhiout(i,j))];
                        S8=['     re-deposited ', num2str(sput_redepos(i,j,p,4)) , ' ; Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        jr=1;
                    end
                    if (jr>npoints)
                        %{
                        errjr=['WARNING, redeposited out of surface area; yr=', num2str(yr), ' jr=', num2str(jr) ' ; setting jr=npoints' ];
                        disp(errjr)
                        S7=['     particle sputtered from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*SpThout_loc(i,j,p)) ' theta =', num2str(180/pi*SpThout(i,j,p)), ' local phi = ' num2str(180/pi*SpPhiout_loc(i,j,p)) ' and phi =', num2str(180/pi*SpPhiout(i,j))];
                        S8=['     re-deposited ', num2str(sput_redepos(i,j,p,4)) , ' ; Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        jr=npoints;
                    end
                    
                    NredepSp(ir,jr)=NredepSp(ir,jr)+1;
                    NSp=NSp+1;
                end
                
            end
            
            %{
            aveSpThout_loc(i,j)=sum(SpPhiout_loc(i,j,p),3)/nEmittSp;
            aveSpThout(i,j)=sum(SpPhiout(i,j,p),3)/nEmittSp;
            aveSpPhiout_loc(i,j)=sum(SpPhiout_loc(i,j,p),3)/nEmittSp;
            aveSpPhiout(i,j)=sum(SpPhiout(i,j,p),3)/nEmittSp;
            %}
            
            aveSpThout_loc(i,j)=aveSpThout_loc_tmp/nEmittSp;
            aveSpThout(i,j)=aveSpThout_tmp/nEmittSp;
            aveSpPhiout_loc(i,j)=aveSpPhiout_loc_tmp/nEmittSp;
            aveSpPhiout(i,j)=aveSpPhiout_tmp/nEmittSp;
            
            
        else %nEmittSp(i,j)==0.0
            
            sput_redepos(i,j,1,1:4)=0.0;
            %S7=['No particles emitted by sputtering from x=', num2str(xi), ', y=', num2str(yj), ' for Y_loc=',num2str(Y_loc(i,j)) , ' and Ncounts=',num2str(Ncounts(i,j)) ];
            %disp(S7)
        end
        
        %%%%%%%%%%%%%%%% 3-Reflection  %%%%%%%%%%%%%%%%
        
        nRpart(i,j)=RN_loc(i,j)*Ncounts(i,j);
        rnR=rand;
       
        if (( nRpart(i,j) - floor(RN_loc(i,j)*Ncounts(i,j)) ) <= rnR )
            nEmittR=floor(nRpart(i,j));
        elseif (( nRpart(i,j) - floor(RN_loc(i,j)*Ncounts(i,j)) ) > rnR )
            nEmittR=floor(nRpart(i,j))+1;
        else
            disp('ERR nEmittSp')
            disp(num2str(nRpart(i,j)))
        end
        
        nRpart(i,j)=nEmittR;
        
        % all reflected particles are emitted with same energy, REout(i,j)
        if(nEmittR>0)
            
            %Stest_nemitt=['for cell in x=', num2str(xi), ' y=', num2str(yj) , ' ; nEmittR=', num2str(nEmittR), ' (RN_loc=',num2str(RN_loc(i,j)) , ' , Ncounts=',num2str(Ncounts(i,j)), ')' ];
            %disp(Stest_nemitt)
            
            %initialize
            RThout(i,j,1:nEmittR)=0.0;
            RThout_loc(i,j,1:nEmittR)=0.0;
            refl_redepos(i,j,1:nEmittR,1:4)=0.0;
            aveRThout_loc_tmp=0.0;
            aveRThout_tmp=0.0;
            
            
            % angle of particles emitted through reflection, wrt local
            % normal to surf:
            RThout_loc(i,j,1:nEmittR) = cosn_distrib(r1, n1, r2, n2, nEmittR, npoints);
            
            %reflected energy and phi angle have a fixed value for each
            %cell, not distributions; define wrt global coordinate system
            REout(i,j)=E0*RE_loc(i,j)/RN_loc(i,j);
            Eout=REout(i,j); %needed?
            
             %specular direction

             %new model:
             % particles moves forward (in x): -pi/2 <= Phi_in <= pi/2
             %if surface_phi also forward -> particle keeps moving in same direction
             %else, apply pure specular reflection:
 
             if (phi>=0 && sg_phi(i,j)>=0 && cos(sg_phi(i,j))>=0) 
                 RPhiout(i,j)=phi;
                 
             elseif(phi<0 && sg_phi(i,j)<0 && cos(sg_phi(i,j))>0) 
                 RPhiout(i,j)=phi;
                 
             else
                 phi_tmp=phi+pi*sign(sg_phi(i,j));
                 RPhiout(i,j)=2*sg_phi(i,j)-phi_tmp;
             end
             
             RPhiout_loc(i,j)=RPhiout(i,j)-sg_phi(i,j);
             
             %{ 
            old model:

            if (sign(phi)==sign(sg_phi(i,j)))
                RPhiout_loc(i,j)=sg_phi(i,j)-phi;
            else
                RPhiout_loc(i,j)=pi+phi-sg_phi(i,j);
            end
          
            RPhiout(i,j) = sg_phi(i,j)-RPhiout_loc(i,j);
            
            if (RPhiout(i,j)>pi)
                RPhiout(i,j)=RPhiout(i,j)-2*pi;
            elseif (RPhiout(i,j)<-pi)
                RPhiout(i,j)=RPhiout(i,j)+2*pi;
            end
            %}
            
            
            %For each reflected particle
            for p=1:nEmittR
                
                %theta angle, wrt global coordinate system
                RThout(i,j,p)=sg_theta(i,j)+RThout_loc(i,j,p);
                %to take averages
                aveRThout_loc_tmp=aveRThout_loc_tmp+RThout_loc(i,j,p);
                aveRThout_tmp=aveRThout_tmp+RThout(i,j,p);
                
                %intersection of reflected particle with surface
                refl_redepos_p=emitted_part_surf_intersec(xi,yj,zk(i,j),RPhiout(i,j),RThout(i,j,p),Ax,fx,Ay,fy,z0,dx,dy);
                refl_redepos(i,j,p,1:4)=refl_redepos_p(1:4); %to store values
                xr=refl_redepos_p(1); %for sorting redep positions
                yr=refl_redepos_p(2);
                zr=refl_redepos_p(3);
                kr=refl_redepos_p(4);
                
                %{
                S7=['     particle reflected from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with theta =', num2str(180/pi*RThout(i,j,p)), ' and phi =', num2str(180/pi*RPhiout(i,j))];
                S8=['     re-deposited: ', num2str(refl_redepos(i,j,p,4)) ,' ;  Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                disp(S7)
                disp(S8)
                fprintf(1, '\n');
                %}
                
                if (kr==0)
                    
                    ir=floor((xr-surfxmin)/di)+1;
                    if (ir<=0)
                        %{
                        errir=['WARNING, redeposited out of surface area; xr=', num2str(xr), ' ir=', num2str(ir) ' ; setting ir=1' ];
                        disp(errir)
                        S7=['     particle reflected from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*RThout_loc(i,j,p)) ' theta =', num2str(180/pi*RThout(i,j,p)), ' local phi = ' num2str(180/pi*RPhiout_loc(i,j)) ' and phi =', num2str(180/pi*RPhiout(i,j))];
                        S8=['     re-deposited: ', num2str(refl_redepos(i,j,p,4)) ,' ;  Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        ir=1;
                    elseif (ir>npoints)
                        %{
                        errir=['WARNING, redeposited out of surface area; xr=', num2str(xr), ' ir=', num2str(ir) ' ; setting ir=npoints' ];
                        disp(errir)
                        S7=['     particle reflected from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*RThout_loc(i,j,p)) ' theta =', num2str(180/pi*RThout(i,j,p)), ' local phi = ' num2str(180/pi*RPhiout_loc(i,j)) ' and phi =', num2str(180/pi*RPhiout(i,j))];
                        S8=['     re-deposited: ', num2str(refl_redepos(i,j,p,4)) ,' ;  Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        ir=npoints;
                    end
                    
                    jr=floor((yr-surfymin)/dj)+1;
                    if (jr<=0)
                        %{
                        errjr=['WARNING, redeposited out of surface area;  yr=', num2str(yr), ' jr=', num2str(jr) ' ; setting jr=1' ];
                        disp(errjr)
                        S7=['     particle reflected from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*RThout_loc(i,j,p)) ' theta =', num2str(180/pi*RThout(i,j,p)), ' local phi = ' num2str(180/pi*RPhiout_loc(i,j)) ' and phi =', num2str(180/pi*RPhiout(i,j))];
                        S8=['     re-deposited: ', num2str(refl_redepos(i,j,p,4)) ,' ;  Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        jr=1;
                    elseif (jr>npoints)
                        %{
                        errjr=['WARNING, redeposited out of surface area; yr=', num2str(yr), ' jr=', num2str(jr) ' ; setting jr=npoints' ];
                        disp(errjr)
                        S7=['     particle reflected from x=', num2str(xi), ', y=', num2str(yj), ', z=', num2str(zk(i,j)), ' with local theta = ' num2str(180/pi*RThout_loc(i,j,p)) ' theta =', num2str(180/pi*RThout(i,j,p)), ' local phi = ' num2str(180/pi*RPhiout_loc(i,j)) ' and phi =', num2str(180/pi*RPhiout(i,j))];
                        S8=['     re-deposited: ', num2str(refl_redepos(i,j,p,4)) ,' ;  Final position: x=', num2str(xr), ', y=', num2str(yr), ', z=', num2str(zr)];
                        disp(S7)
                        disp(S8)
                        %}
                        jr=npoints;
                    end
                    
                    NredepR(ir,jr)=NredepR(ir,jr)+1;
                    NR=NR+1;
                end
                
            end
            
            aveRThout_loc(i,j)=aveRThout_loc_tmp/nEmittR;
            aveRThout(i,j)=aveRThout_tmp/nEmittR;
            
            %aveRThout_loc(i,j)=sum(RThout_loc(i,j,p),3)/nEmittR;
            %aveRThout(i,j)=sum(RThout(i,j,p),3)/nEmittR;
            
            
        else %nEmittR(i,j)==0.0
            
            refl_redepos(i,j,1,1:4)=0.0;
            %S7=['No particles emitted by reflection from x=', num2str(xi), ', y=', num2str(yj), ' for Y_loc=',num2str(RN_loc(i,j)) , ' and Ncounts=',num2str(Ncounts(i,j)) ];
            %disp(S7)
        end
        %{
        if (nEmittSp>0 || nEmittR>0)
            Sredep=['    Number of particles redeposited:  from sputtering ',num2str(NSp), ', from reflection ' , num2str(NR) ];
            disp(Sredep)
            fprintf(1, '\n');
        end
        %}
        
        
    end
    
    if (mod(10*i/npoints,1)==0)
        S4=['   ...',num2str(100*i/npoints),'% done'];
        disp(S4)
    end
end

fprintf(1, '\n');
disp('---- particle emission done -----')
fprintf(1, '\n');

save(filename,'NredepSp', '-append');
save(filename,'NredepR', '-append');

save(filename,'nSppart', '-append');
save(filename,'SpEout','-append');
save(filename,'SpThout','-append');
save(filename,'SpThout_loc','-append');
save(filename,'aveSpThout','-append');
save(filename,'aveSpThout_loc','-append');
save(filename,'SpPhiout', '-append');
save(filename,'SpPhiout_loc', '-append');
save(filename,'aveSpPhiout', '-append');
save(filename,'aveSpPhiout_loc', '-append');
save(filename,'sput_redepos', '-append');

save(filename,'nRpart', '-append');
save(filename,'REout','-append');
save(filename,'RThout','-append');
save(filename,'RThout_loc','-append');
save(filename,'aveRThout','-append');
save(filename,'aveRThout_loc','-append');
save(filename,'RPhiout', '-append');
save(filename,'RPhiout_loc', '-append');
save(filename,'refl_redepos', '-append')

