using Rotations
using LinearAlgebra


function trans(x)
    
    alpha=x[1];  a=x[2]; d=x[3]; t=x[4]

   T=[ cos(t)                       -sin(t)          0              a;
    sin(t)*cos(alpha)         cos(t)*cos(alpha) -sin(alpha)  -sin(alpha)*d;
    sin(t)*sin(alpha)         cos(t)*sin(alpha)  cos(alpha)   cos(alpha)*d;
    0                          0                         0          1     ]
    return T
    
end

function J(t)
    t1= t[1] ; t2= t[2];
    t3= t[3]; 
    t4= t[4];
    t5= t[5]; 
    t6= t[6];
    j=zeros(3,6)
        j[1,1]= 0.227*sin(t1)*sin(t5)*cos(t2 + t3 + t4) - 0.095*sin(t1)*sin(t2 + t3 + t4) + 0.425*sin(t1)*cos(t2) + 0.392*sin(t1)*cos(t2 + t3) + 0.227*cos(t1)*cos(t5) + 0.109*cos(t1)       
        j[1,2]= 0.425*sin(t2)*cos(t1) + 0.227*sin(t5)*sin(t2 + t3 + t4)*cos(t1) + 0.392*sin(t2 + t3)*cos(t1) + 0.095*cos(t1)*cos(t2 + t3 + t4)       
        j[1,3]= 0.227*sin(t5)*sin(t2 + t3 + t4)*cos(t1) + 0.392*sin(t2 + t3)*cos(t1) + 0.095*cos(t1)*cos(t2 + t3 + t4)       
        j[1,4]= 0.227*sin(t5)*sin(t2 + t3 + t4)*cos(t1) + 0.095*cos(t1)*cos(t2 + t3 + t4)       
        j[1,5]= -0.227*sin(t1)*sin(t5) - 0.227*cos(t1)*cos(t5)*cos(t2 + t3 + t4)       
        j[1,6]= 0       
        
        j[2,1]= 0.227*sin(t1)*cos(t5) + 0.109*sin(t1) - 0.227*sin(t5)*cos(t1)*cos(t2 + t3 + t4) + 0.095*sin(t2 + t3 + t4)*cos(t1) - 0.425*cos(t1)*cos(t2) - 0.392*cos(t1)*cos(t2 + t3)       
        j[2,2]= 0.425*sin(t1)*sin(t2) + 0.227*sin(t1)*sin(t5)*sin(t2 + t3 + t4) + 0.392*sin(t1)*sin(t2 + t3) + 0.095*sin(t1)*cos(t2 + t3 + t4)       
        j[2,3]= 0.227*sin(t1)*sin(t5)*sin(t2 + t3 + t4) + 0.392*sin(t1)*sin(t2 + t3) + 0.095*sin(t1)*cos(t2 + t3 + t4)       
        j[2,4]= 0.227*sin(t1)*sin(t5)*sin(t2 + t3 + t4) + 0.095*sin(t1)*cos(t2 + t3 + t4)       
        j[2,5]= -0.227*sin(t1)*cos(t5)*cos(t2 + t3 + t4) + 0.227*sin(t5)*cos(t1)       
        j[2,6]= 0       
        
        j[3,1]= 0       
        j[3,2]= -0.227*sin(t5)*cos(t2 + t3 + t4) + 0.095*sin(t2 + t3 + t4) - 0.425*cos(t2) - 0.392*cos(t2 + t3)       
        j[3,3]= -0.227*sin(t5)*cos(t2 + t3 + t4) + 0.095*sin(t2 + t3 + t4) - 0.392*cos(t2 + t3)       
        j[3,4]= -0.227*sin(t5)*cos(t2 + t3 + t4) + 0.095*sin(t2 + t3 + t4)       
        j[3,5]= -0.227*sin(t2 + t3 + t4)*cos(t5)       
        j[3,6]= 0     


        return j
end

function CMP(e)
    pos=[]
    for i in range(1,size(e)[1])
        te=e[i,1:3];  eff_ang=e[i,4:6]';
        re = round.(RotXYZ(eff_ang[1],eff_ang[2],eff_ang[3]),digits=3);
        Xwe=(vcat(hcat(re,te),[0 0 0 1]));
            function CM(t)
                t1= t[1] ; t2= t[2];
                t3= t[3]; 
                t4= t[4];
                t5= t[5]; 
                t6= t[6]; 
                Xbe= Matrix(I,4,4) # e's coordinate in b's frame
                a2=-.425 ; d3=-.120 ; a3=-0.392; d4=0.093; # all are in mm
                dh=[0      -.1395    .482      t1;
                    pi/2     0        0.136      t2;
                    0        a2       d3       t3;
                    0        a3       d4       t4;
                    pi/2     0       .095      t5;
                -pi/2     0       .227      t6;] 
                
                for i in 1:(size(dh)[1])
                    # display(trans(dh[i,:]))
                    Xbe= Xbe*trans(dh[i,:])
                end     
                return Xwe*inv(Xbe), Xbe[1:3,1:3]  ;   
            end

            T=[]
            for t1 in range(0,360,15)
                for t2 in range(0,180,15)
                    for t3 in range(0,180,15)
                        for t4 in range(0,180,15)
                            for t5 in range(0,360,15)
                                for t6 in range(0,360,15)
                                    t=[t1 t2 t3 t4 t5 t6];
                                    p= CM(t);
                                    c=Rotations.params(RotXYZ(p[2]))[1:2,1];                        
                                    for i in range(1,2)
                                        if c[i]<0
                                            c[i]=2*pi+c[i];
                                        end
                                    end                            
                                    if  0<abs(c[1]-eff_ang[1])<0.05 && 0<abs(c[2]-eff_ang[2])<0.05 && 0<=p[1][15]< 0.0025
                                        l=J(t);
                                        u=sqrt(abs(det(l*l')));
                                        v=(p[1],u,t);
                                        T=vcat(T,v);
                                    end                                                 
                                end
                            end
                        end
                    end
                end
            end
            
        pos=vcat(pos,Tuple(T));
    end
    return pos

end
e2=[        0.0       0.0          0.2          1.3 0.4 0.8;
            2.51327   0.951057     0.407295     1.5 0.3 0.7;
            5.02655   0.587785     0.742705     1.0 0.5 0.8;
            7.53982  -0.587785     0.742705     2.3 1.0 1.2;
            10.0531   -0.951057     0.407295    0.56 1.5 1.0;
            12.5664   -2.44929e-16  0.2         1.2 0.2 1.2;
];

cmpalg_1=CMP(e2)