function T = rotation_material_matrix( alpha,beta,gamma )
% get the rotation matrix for material matrix
a=deg2rad(alpha);
b=deg2rad(beta);
r=deg2rad(gamma);

% Z-Y-X Rotate
Rx=[1,0,0;0, cos(a), sin(a);0,-sin(a),cos(a)];
Ry=[cos(b),0, -sin(b);0,1,0; sin(b),0,cos(b)];
Rz=[cos(r),sin(r),0;-sin(r),cos(r),0;0,0,1];
a =Rx*Ry*Rz;
T(1:6,1:6) = 0;
for i=1:1:3
    for j=1:1:3
        if i==j
            alpha = j; 
        else
            alpha = 9-i-j; 
        end
        for p=1:1:3
            for q=1:1:3
                if p==q 
                    beta = p; 
                else
                    beta = 9-p-q; 
                end
                T(alpha,beta) = 0;
                if alpha<=3 && beta<= 3 
                    T(alpha,beta)=a(i,p)*a(i,p); 
                end
                if alpha> 3 && beta<= 3
                    T(alpha,beta)=a(i,p)*a(j,p); 
                end
                if alpha<=3 && beta>3
                    T(alpha,beta)=a(i,q)*a(i,p)+a(i,p)*a(i,q);
                end
                if alpha>3 && beta>3
                    T(alpha,beta)=a(i,p)*a(j,q)+a(i,q)*a(j,p);
                end
            end
        end
    end
end
%% Finite Element Analysis of Composite Materials Using Abaqus P14
% Derivation of the transformation matrix [T]
% clear all;
% syms T alpha R
% syms a a11 a12 a13 a21 a22 a23 a31 a32 a33
% a = [a11,a12,a13;
% a21,a22,a23;
% a31,a32,a33];
% T(1:6,1:6) = 0;
% for i=1:1:3
% for j=1:1:3
% if i==j; alpha = j; else alpha = 9-i-j; end
% for p=1:1:3
% for q=1:1:3
% if p==q beta = p; else beta = 9-p-q; end
% T(alpha,beta) = 0;
% if alpha<=3 & beta<= 3; T(alpha,beta)=a(i,p)*a(i,p); end
% if alpha> 3 & beta<= 3; T(alpha,beta)=a(i,p)*a(j,p); end
% if alpha<=3 & beta>3; T(alpha,beta)=a(i,q)*a(i,p)+a(i,p)*a(i,q);end
% if alpha>3 & beta>3; T(alpha,beta)=a(i,p)*a(j,q)+a(i,q)*a(j,p);end
% end
% end
% end
% end
% T
% R = eye(6,6); R(4,4)=2; R(5,5)=2; R(6,6)=2; % Reuter matrix
% Tbar = R*T*R^(-1)

%% PhD Thesis of Di shengjie; I think it may be wrong
% l1=sin(a)*cos(b)*cos(r)+cos(a)*sin(r);
% m1=cos(a)*cos(b)*cos(r)-sin(a)*sin(r);
% n1=-sin(b)*cos(r);
% 
% l2=sin(a)*cos(b)*sin(r)-cos(a)*cos(r);
% m2=cos(a)*cos(b)*sin(r)+sin(a)*cos(r);
% n2=-sin(b)*sin(r);
% 
% l3=sin(a)*sin(b);
% m3=cos(a)*sin(b);
% n3=cos(b);
% 
% T=[l1^2,l2^2,l3^3,2*l2*l3,2*l1*l3,2*l1*l2;...
%     m1^2,m2^2,m3^3,2*m2*m3,2*m1*m3,2*m1*m2;...
%     n1^2,n2^2,n3^3,2*n2*n3,2*n1*n3,2*n1*n2;...
%     m1*n1,m2*n2,m3*n3,m2*n3+m3*n2,m1*n3+m3*n1,m1*n2+m1*n1;...
%     n1*l1,n2*l2,n3*l3,n2*l3+n3*l2,n1*l3+n3*l1,n1*l2+n1*l1;...
%     m1*l1,m2*l2,m3*l3,m2*l3+m3*l2,m1*l3+m3*l1,m1*l2+m1*l1];

end

