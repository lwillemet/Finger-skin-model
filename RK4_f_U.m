function dU = f_U(tc,Uc,Fc)
global km kt kn Kn B posi Ndx nu

pos_x = Uc(Ndx+2:2*(Ndx+1))+posi(Ndx+2:2*(Ndx+1));
pos_z = Uc(1:Ndx+1)+posi(1:Ndx+1);

part_int = [-ones(1,floor(Ndx/2)),0,ones(1,floor(Ndx/2))];
part_ext = [-ones(1,floor(Ndx/2)),ones(1,floor(Ndx/2))];

    %%Stiffness matrix (spring)
% % Km: Springs dependencies of the membrane
angleA = atan(diff(pos_z(1:Ndx+1))./diff(pos_x(1:Ndx+1))); %% angle between 2 elements and the horizontal spring
Acos = abs(cos(angleA(1:Ndx)))';
Asin = abs(sin(angleA(1:Ndx)))';

% Boundary conditions
Acos(1,1) = 0;
Asin(1,1) = 0;
Acos(1,Ndx+1) = 0;
Asin(1,Ndx+1) = 0;

d0 = -[0,Asin(1:Ndx)+Asin(2:Ndx+1),0,Acos(1:Ndx)+Acos(2:Ndx+1)]';
d0sub = [0,Asin(2:Ndx),0,0,Acos(2:Ndx),0]';        %% sous-diagonale
d0sup = [0,0,Asin(2:Ndx),0,0,Acos(2:Ndx)]';        %% sur-diagonale
d12 = -[0,[0,part_ext].*Acos(1:Ndx)+[part_ext,0].*Acos(2:Ndx+1)]';
d12sub = [0,part_ext.*Acos(2:Ndx),0]';
d12sup = [0,0,part_ext.*Acos(2:Ndx)]';
d21 = -[0,[0,part_ext].*Asin(1:Ndx)+[part_ext,0].*Asin(2:Ndx+1)]';
d21sub = [0,part_ext.*Asin(2:Ndx),0]';
d21sup = [0,0,part_ext.*Asin(2:Ndx)]';
spKm = (spdiags([d0sub d0 d0sup],-1:1,length(d0),length(d0)))+...
    0*(spdiags([zeros(Ndx+1,3);[d12sub d12 d12sup]],Ndx:Ndx+2,length(d0),length(d0)))+...
    0*(spdiags([[d21sub d21 d21sup];zeros(Ndx+1,3)],-Ndx-2:-Ndx,length(d0),length(d0)));

% % Kt: Springs dependencies of the subcutaneous tissues
angleB = atan(abs(pos_z(1)-pos_z(2:Ndx+1))./abs(pos_x(1)-pos_x(2:Ndx+1))); %% angle between two elements linked to the bone
Bcos = abs(cos(angleB(1:Ndx)))';
Bsin = abs(sin(angleB(1:Ndx)))';

v11 = [-sum(Bsin),Bsin];
v22 = [-sum(Bcos),Bcos];
v12 = [sum(part_int.*Bcos(1:Ndx)),-part_int.*Bcos(1:Ndx)];
v21 = [sum(part_int.*Bsin(1:Ndx)),-part_int.*Bsin(1:Ndx)];

spKt = (sparse(1,1:Ndx+1,v11,2*(Ndx+1),2*(Ndx+1))+...
    sparse(2:Ndx+1,1,Bsin,2*(Ndx+1),2*(Ndx+1))+...
    sparse(Ndx+2,Ndx+2:2*(Ndx+1),v22,2*(Ndx+1),2*(Ndx+1))+...
    sparse(Ndx+3:2*(Ndx+1),Ndx+2,Bcos,2*(Ndx+1),2*(Ndx+1))+...
    spdiags([0,-Bsin,0,-Bcos]',0,2*(Ndx+1),2*(Ndx+1)))+...
    nu*(sparse(1,Ndx+2:2*(Ndx+1),v12,2*(Ndx+1),2*(Ndx+1))+...
    sparse(2:Ndx+1,Ndx+2,-part_int.*Bcos,2*(Ndx+1),2*(Ndx+1))+...
    spdiags([zeros(1,Ndx+1),0,part_int.*Bcos]',Ndx+1,2*(Ndx+1),2*(Ndx+1)))+...
    nu*(sparse(Ndx+2,1:Ndx+1,v21,2*(Ndx+1),2*(Ndx+1))+...
    sparse(Ndx+3:2*(Ndx+1),1,-part_int.*Bsin,2*(Ndx+1),2*(Ndx+1))+...
    spdiags([0,part_int.*Bsin]',-Ndx-1,2*(Ndx+1),2*(Ndx+1)));

% % Sum of springs dependancies matrix
spK = km.*spKm+kt.*spKt+kn.*Kn;
 
dU = -inv(B)*spK*Uc-inv(B)*Fc;

end



