using DataFrames
using Pandas
using PyCall
using PyPlot
plt=pyimport("matplotlib.pyplot")
np=pyimport("numpy")
df = DataFrames.DataFrame(read_json("cmpalg_1.json"));
sz=size(df,2)
global P1=[]; global P2=[]; global  P3=[]; global P4=[]; 
global P5=[]; global P6=[]; 
e2=[        0.0       0.0          0.2         ;
            2.51327   0.951057     0.407295    ;
            5.02655   0.587785     0.742705     ;
            7.53982  -0.587785     0.742705    ;
            10.0531   -0.951057     0.407295    ; 
            12.5664   -2.44929e-16  0.2       ;
];
function srt(a)
    if size(a,1)<15
        s = size(a,1)
    else
        s = Int32(round(0.2*size(a,1)))
    end
    return s
end

for i in range(1,sz)
    if df[1,i]==nothing
        break
    end
    global P1 = vcat(P1,Tuple(df[1,i]))
end
for i in range(1,sz)
    if df[2,i]==nothing
        break
    end
    global P2 = vcat(P2,Tuple(df[2,i]))
end
for i in range(1,sz)
    if df[3,i]==nothing
        break
    end
    global P3 = vcat(P3,Tuple(df[3,i]))
end
for i in range(1,sz)
    if df[4,i]==nothing
        break
    end
    global P4 = vcat(P4,Tuple(df[4,i]))
end
for i in range(1,sz)
    if df[5,i]==nothing
        break
    end
    global P5 = vcat(P5,Tuple(df[5,i]))
end
for i in range(1,sz)
    if df[6,i]==nothing
        break
    end
    global P6 = vcat(P6,Tuple(df[6,i]))
end
pb1=sort!(P1,by=x->x[2],rev=true);
pb2=sort!(P2,by=x->x[2],rev=true);
pb3=sort!(P3,by=x->x[2],rev=true);
pb4=sort!(P4,by=x->x[2],rev=true);
pb5=sort!(P5,by=x->x[2],rev=true);
pb6=sort!(P6,by=x->x[2],rev=true);
fig=plt.figure(figsize=(10,10))
ax = plt.axes(projection ="3d")
pbu1=pb1[1:srt(pb1)];pbu2=pb2[1:srt(pb2)]
pbu3=pb3[1:srt(pb3)]
pbu4=pb4[1:srt(pb4)]
pbu5=pb5[1:srt(pb5)]
pbu6=pb6[1:srt(pb6)]
function plotp(p)
    mat=Array{Float64}(undef, 3, 0)
    for i in range(1,size(p,1))
        mat = hcat(mat,p[i][1][1:3,4])
    end
    ax.scatter3D(mat[1,:],mat[2,:], zeros(1,size(mat,2)),label="BaseCoordinates")
    xlim(0,15); ylim(-5,5); zlim(0,2);
    return mat
end
plotp(pb1); plotp(pb2); plotp(pb3); plotp(pb4); plotp(pb5); plotp(pb6)
ax.scatter3D(e2[:,1],e2[:,2],e2[:,3],s=155,cmap="rainbow")
plt.show()



fig2=plt.figure(figsize=(50,50))
ax = plt.axes(projection ="3d")
q1=plotp(pbu1); q2=plotp(pbu2); q3=plotp(pbu3); q4=plotp(pbu4);q5=plotp(pbu5); q6=plotp(pbu6)
np1 = np.sum(q1,axis=1)/size(q1,2)
np2 = np.sum(q2,axis=1)/size(q2,2)
np3 = np.sum(q3,axis=1)/size(q3,2)
np4 = np.sum(q4,axis=1)/size(q4,2)
np5 = np.sum(q5,axis=1)/size(q5,2)
np6 = np.sum(q6,axis=1)/size(q6,2)

npm= [np1 np2 np3 np4 np5 np6]
ax.scatter3D(e2[:,1],e2[:,2],e2[:,3],s=155,cmap="rainbow",label="Eff Trajectory")
ax.scatter3D(npm[1,:],npm[2,:], zeros(1,size(npm,2)),color="k",s=100, label="Optimum_Points")
ax.legend()
plt.show()

