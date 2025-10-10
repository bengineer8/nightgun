#version 460 core

#define positionsSize 1048576

const float tpi=6.283185307;//two pi
const int shaderDataLength=8*128+8*4096+4*8192;
const int sizeOfpow = 14;
const float pr=1.0/9.0, prc=cos(pr);
const float wallthickness=1/18.0;
const float wallcornershad=1.0/9;
const float coswallthickness=cos(wallthickness);
const float coswallcornershad=cos(wallcornershad);
const float inf = 1.0/0;
//const vec4 bgcolor=vec4(0,.125,0.498,1);
const vec4 bgcolor=vec4(0,0,0,1);
const int cellCount = 211;
out vec4 FragColor;
uniform int res=800;
uniform float T=1.0,pz,pcx,pcy,pcz,zoom=1.0;
uniform vec3 pl,pl2,camRef,duppl,duppl2;
uniform int pw;//player world
uniform int duppw;
uniform float mirrorPlayer=1;//-1 for mirror player

layout (std430, binding=3) buffer ssbo1
{
  float positions[positionsSize];
};
uniform float iaunits[16];

float Mod(float a, float b){
    a = mod(a,b);
    if(a < 0) a = b + a;
    if(a == b) a = 0;
    return a;
}

float floorlog2(float f){
    uint I = floatBitsToUint(f);
    I = I>>23;
    I = I&0xff;
    return float(I) - 127;
}

float exponentialFloorb2(float f){//largest power of 2 less than f, awaiting a better name
    uint I = floatBitsToUint(f);
    I = I&(0xff800000);
    return uintBitsToFloat(I);
}

mat3 h2invert(mat3 a){//if a is a combination of rotations and translations in h2, b is the inverse
    mat3 b=transpose(a);
    b[2][0]=-b[2][0];
    b[2][1]=-b[2][1];
    b[0][2]=-b[0][2];
    b[1][2]=-b[1][2];
    return b;
}

float innerdot(vec3 a, vec3 b){
    return dot(a,b*vec3(1,1,-1));
}

vec3 lcross(vec3 a, vec3 b){
    return cross(a,b)*vec3(1,1,-1);
}

float s2disrank(float z, float x){//gives a quick metric for the distance in a specific direction from (0,0,1)
    z+=1;
    if(x<0) z=-z;
    return 2-z;
}

float distance2(vec2 a, vec2 b){
    vec2 c = a - b;
    return dot(c,c);
}

float distance2(vec3 a, vec3 b){
    return dot(a-b,a-b);
}

float arctan(float x, float z){//makes sense in context, with me wanting the distance from (0,0,1) in a specific direction
    x=atan(x,z);
    if(x<0) x+=tpi;
    return x;
}

float smallArccos(float x){//small angle approximation for acos
    return sqrt(2 - 2*x);
}

float safe_sqrt(float x){//limited sqrt
    if(x<0) return 0;
    return sqrt(x);
}

mat3 s2matto(vec3 p){//returns a matrix that rotates a unit sphere so p is at (0,0,1)
    float r=sqrt(p.x*p.x+p.y*p.y), c=p.x/r, s=p.y/r;
    if(r==0){
        c=1;s=0;
    }
    return mat3(p.z*c,-s,p.x, p.z*s,c,p.y, -r,0,p.z);
}

mat3 s2matto(int i){
    float x = positions[i],
    y = positions[i + 1],
    z = positions[i + 2],
    ns = positions[i + 8],
    c = positions[i + 9],
    nr = positions[i + 10],
    c2 = positions[i + 11],
    s2 = positions[i + 12];
    return mat3(c2,-s2,0, s2,c2,0, 0,0,1)*mat3(z*c,ns,x, -z*ns,c,y, nr,0,z);
    //return mat3(1,0,0, 0,1,0, 0,0,1);
}

bool e2pbc(vec3 pi, vec3 pf, int w){//e2 portal between check, checks if there is a portal between 2 points in world w
    bool yes=false;//coding is my passion
    mat3 rot;
    float d;
    {
        vec3 v=(pf-pi);
        d=length(v);
        v/=d;
        rot=mat3(v.x,-v.y,0,v.y,v.x,0,0,0,1);
    }
    int I=int(positions[w]);
    while(I<int(positions[w+1])&&!yes){
        vec3 p1=vec3(positions[I],positions[I+1],0), p2=vec3(positions[I+3],positions[I+4],0);
        float r=distance(p1,p2);
        int type=int(positions[I+7])&3;
        float al=positions[I+6];//angle limit, not an angle.
        p1=rot*(p1-pi);p2=rot*(p2-pi);//this is the reverse order it is done in the main function
        if(type==0 && abs(p1.y)<r && (p2.y*p2.y<=al*r*r || al==4)){
            float o=sqrt(r*r-p1.y*p1.y), xt=p1.x-o;
            if(((al==4) || ((xt-p2.x)*(xt-p2.x)+p2.y*p2.y<al*r*r)) && xt>0 && xt<d) {
                yes=true;
            }
            xt=p1.x+o;
            if(((al==4) || ((xt-p2.x)*(xt-p2.x)+p2.y*p2.y<al*r*r)) && xt>0 && xt<d) {
                yes=true;
            }
        }
        if(type==1 && p1.y*p2.y<=0 && (p1.x>0||p2.x>0)){
            float xt=(p1.x-p2.x)*p2.y/(p2.y-p1.y)+p2.x;
            if(xt>0 && xt<d){
                yes=true;
            }
        }
        I += sizeOfpow;
    }
    return yes;
}//a lot of this is copying code from earlier. Might generalize later and spin out the reused code as its own function

bool s2pbc(vec3 pi, vec3 pf, int w){//s2 portal between check
    bool yes=false;//:P
    mat3 rot=s2matto(pi);
    pf=rot*pf;
    float c=sqrt(pf.x*pf.x+pf.y*pf.y), s;
    float disrank=s2disrank(pf.z,c);
    c=1/c;
    s=pf.y*c;
    c*=pf.x;
    rot=mat3(c,-s,0, s,c,0, 0,0,1)*rot;
    int I=int(positions[w]);
    while(I<int(positions[w+1]) && !yes){
        vec3 p1=vec3(positions[I],positions[I+1],positions[I+2]), p2=vec3(positions[I+3],positions[I+4],positions[I+5]);
        float rc=dot(p1,p2);
        float rs=1-rc*rc;//actually rs^2, just reusing space
        p1=rot*p1;
        p2=rot*p2;
        if(p1.y*p1.y<rs){
            vec3 cipp, cipm;//closest intersection point plus and minus
            float s=sqrt(rs-p1.y*p1.y), D=1/(1-p1.y*p1.y);
            cipp=vec3(p1.x*rc+p1.z*s,0,p1.z*rc-p1.x*s)*D;
            cipm=vec3(p1.x*rc-p1.z*s,0,p1.z*rc+p1.x*s)*D;
            float al=positions[I+6];//angle limit, (cos(a)-1), decreases with angle
            float ala=1+rs*al;
            if( (al==-2 || dot(cipp,p2)>ala ) && s2disrank(cipp.z,cipp.x)<disrank){//either not the source portal or the other point is closer
                yes=true;
            }
            if( (al==-2 || dot(cipm,p2)>ala ) && s2disrank(cipm.z,cipm.x)<disrank){//either not the source portal or the other point is closer
                yes=true;
            }
        }
        I += sizeOfpow;
    }
    return yes;
}


mat3 h2matto(vec3 p){
    float r = sqrt(p.x*p.x + p.y*p.y), c = p.x/r, s = p.y/r;
    if(r == 0){
        c = 1;
        s = 0;
    }
    return mat3(p.z*c,-s,-p.x, p.z*s,c,-p.y, -r,0,p.z);
}

void main(){
    int wallIndexes[64];
    int current_world=pw;
    vec2 screenPos=vec2(2.0*gl_FragCoord.x/float(res)-1.0,2.0*gl_FragCoord.y/float(res)-1.0);
    float t,I,od,mirror=mirrorPlayer;//make mirror -1 for mirror
    if(length(screenPos)<1.0){
        screenPos.y*=mirrorPlayer;
        int type=0,si=-1;//source index
        float d=length(screenPos);
        //float d=length(screenPos), dat=33554432;//impossible hard coded value as a signal that the source is the player, 2^25
        screenPos=screenPos/d;
        d=d*tpi*zoom;
        od=d;
        vec3 sp1 = pl, sp2=camRef;//source p1 and p2
        mat3 lm = mat3(screenPos[0],-screenPos[1],0,screenPos[1],screenPos[0],0,0,0,1),am;//loaction matrix and angle matrix (unused when player is the source)
        mat3 rotstart;
        if(iaunits[current_world] > 0){
                rotstart = s2matto(sp1);
                sp2 = rotstart*sp2;
                t = sqrt(sp2.x*sp2.x + sp2.y*sp2.y);
                sp2 /= t;
                rotstart = mat3(sp2.x,-sp2.y,0, sp2.y,sp2.x,0, 0,0,1)*rotstart;
        }
        vec3 endOfRay;
        float iaunit;
        bool hitwall = false;
        int iterations = 0;
        while(d>0.0&&iterations<9){
            iaunit = iaunits[current_world];
            iterations++;//saftey
            if(iaunit == 0 && iterations < 9){//E2
                vec3 off=-sp1;//offset
                sp2=sp2+off;
                t = length(sp2);
                sp2 /= t;
                mat3 rot=mat3(sp2[0],-sp2[1],0,sp2[1],sp2[0],0,0,0,1);
                if(type==0) rot=lm*rot;
                off=rot*off;//keeping track of where (0,0) ends up for later.
                if(si>-1){//accounts for the angle the portal was hit at. Disabled when the source is the player
                    if(type==0) off.x-=t;
                    if(type==1) {
                        off.x-=lm[0][1]+t/2;
                        am=mat3(0,1,0, -1,0,0, 0,0,1)*am;
                    }
                    off=am*off;
                    rot=am*rot;
                }
                mat3 roti=transpose(rot);//quick inverse
                //to convert from global to local coords: rot*(p)+off, the reverse is roti*(p-off)
                int i=int(positions[current_world]),cii=-1;
                vec3 cip1,cip2;
                float cid=d;
                wallIndexes[0]=0;
                while(i<int(positions[current_world+1]) && iterations<9){
                    vec3 p1=vec3(positions[i],positions[i+1],0), p2=vec3(positions[i+3],positions[i+4],0);
                    t=distance(p1,p2);
                    type = int(positions[i+7])&7;
                    if(type > 3){
                        type -= 4;
                        wallIndexes[0]++;
                        wallIndexes[wallIndexes[0]]=i;
                    }
                    float al=positions[i+6];//angle limit, not an angle.
                    p1=rot*(p1)+off;p2=rot*(p2)+off;
                    if(type == 0 && abs(p1.y)<t&&(p2.y*p2.y<=al*t*t || al==4)){
                        float o=sqrt(t*t-p1.y*p1.y), xt=p1.x-o, xt1=xt;
                        if(((al==4) || ((xt-p2.x)*(xt-p2.x)+p2.y*p2.y<al*t*t)) && xt>0 && xt<cid && (i!=si)) {//1st collision is not allowed if portal is the source
                            cip1=p1;cip2=p2;cii=i;cid=xt;
                        }
                        xt=p1.x+o;
                        if(((al==4) || ((xt-p2.x)*(xt-p2.x)+p2.y*p2.y<al*t*t)) && xt>0 && xt<cid && (i!=si||abs(xt1)<xt)) {//2nd clollision is only allowed with oneself if the 1st collsion would be the starting point
                            cip1=p1;cip2=p2;cii=i;cid=xt;
                        }
                    }
                    if(type == 1 && i!=si && p1.y*p2.y<=0 && (p1.x>0 || p2.x>0) ){
                            float xt=(p1.x-p2.x)*p2.y/(p2.y-p1.y)+p2.x;
                            if(xt>0 && xt<cid){
                                cip1=p1;cip2=p2;cii=i;cid=xt;
                            }
                        }
                    i += sizeOfpow;
                }
                if(cii>-1) {
                    int dat=int(positions[cii+7]);
                    if( (dat&4)>0 ) {
                        FragColor=vec4(0,0,0,1);
                        hitwall = true;
                        d=0;
                    } else {
                        float side=1;
                        if(((dat>>22)&1)==1) side=-1;
                        mirror=1;
                        if(((dat>>23)&1)==1) mirror=-1;
                        type=dat&3;
                        current_world=(dat>>3)&511;
                        int di=int(positions[current_world])+sizeOfpow*(dat>>12)&1023;
                        sp1=vec3(positions[di],positions[di+1],positions[di+2]);
                        sp2=vec3(positions[di+3],positions[di+4],positions[di+5]);
                        si=di;
                        t=distance(cip1,cip2);
                        cip1.x-=cid;
                        if(type == 0){
                            cip2=(cip2-vec3(cid,0,0)-cip1)/t;
                            cip1=-cip1/t;//just reusing the variable
                            am=mat3(cip1.x,mirror*cip1.y,0,-mirror*cip1.y,cip1.x,0,0,0,side*1)*side;
                            cip1=mat3(cip2.x,-cip2.y,0,cip2.y,cip2.x,0,0,0,1)*cip1;
                            lm=mat3(cip1.x,mirror*-cip1.y*side,0,mirror*cip1.y*side,cip1.x,0,0,0,1);
                        }
                        else{
                            float t2=length(cip1);
                            cip1.x*=mirror;
                            cip1*=side;
                            am=mat3(-cip1.y,cip1.x,0, -cip1.x,-cip1.y,0, 0,0,t2)/t2;
                            //lm[1][0]=(t2-t/2)*mirror*side;
                            float sL = -(t2-t/2)*mirror*side;
                            lm = mat3(1 - sL*sL/2,-sL,-sL*sL/2, sL,1,sL,  sL*sL/2,sL,1 + sL*sL/2);//NOTE: this may break something!
                        }
                        d-=cid;
                        //FragColor=vec4(1-(d-cid)*abs(am[0][0]),0,0,1);d=0;

                        //d=0.0;//debug view showing their relative coordinates
                        //FragColor = vec4(lm[0][0],0.0*float(cii)/float(I),lm[1][0],1.0);
                    }
                } else {//didn't hit anything, local to global coords now
                    endOfRay = roti*(vec3(d,0,0) - off);
                    /*
                    //
                    vec3 gl=roti*(vec3(d,0,0)-off);
                    vec2 cp, ccp, RL;
                    i=1;
                    int iwi=-1;
                    vec2 cp1, cp2;
                    float crl;
                    float ambientoc=0, cd=2*wallthickness*wallthickness;
                    while(i<wallIndexes[0]+1){
                        int wi=wallIndexes[i];
                        vec2 p1=vec2(positions[wi],positions[wi+1]);
                        vec2 p2=vec2(positions[wi+3],positions[wi+4])-p1;
                        vec2 glc=vec2(gl.x,gl.y)-p1;
                        float rl=dot(p2,p2);//radius/ length squared
                        type=int(positions[wi+7])&3;
                        if(type==0){
                            cp=glc;
                            float al=positions[wi+6], ala=al*rl;
                            t=sqrt(rl/(dot(glc,glc)));
                            cp*=t;
                            RL=cp-p2;
                            if(al<4 && dot(RL,RL)>ala){
                                float c=1-al*0.5, s=sqrt(1-c*c);
                                if(glc.y*p2.x-glc.x*p2.y<0) s=-s;
                                cp=mat2(c,s, -s,c)*p2;
                            }
                        } else {
                            t=min(max(0,dot(glc,p2)),rl);
                            cp=p2*t/rl;
                        }
                        RL=cp-glc;
                        t=dot(RL,RL);
                        if(t<pr*pr){
                            if(t<cd*.99){
                                cd=t;
                                if(t<wallthickness*wallthickness) {
                                    iwi=wi;
                                    cp1=p1;
                                    cp2=p2;
                                    crl=rl;
                                    ccp=cp;
                                }
                                else ambientoc=t;
                            }
                        }
                        i+=1;
                    }*/
                    ////////////////////
                    /*if(iwi<0){
                        if(ambientoc>0) ambientoc=(sqrt(ambientoc)/pr)/2+0.5;
                        else ambientoc=1;
                        float zx=gl.x,zy=gl.y;
                        int N=60,n=N;
                        while((n>0)&&(zx*zx+zy*zy<4.0)){//mset as a hello world I nver removed. It looks cool.
                            t=zx*zx-zy*zy+gl.x;
                            zy=2.0*zx*zy+gl.y;
                            zx=t;
                            n--;
                        }
                        I=float(n)/float(N);
                        if(distance(gl,pl)<pr && current_world==pw && !e2pbc(gl,pl,pw)){
                            float grad=1-distance(gl,pl)/pr;
                            if(distance(gl,pl2)<pr) FragColor=vec4(1,grad,1,1);
                            else FragColor=vec4(0,grad,1,1);
                        } else if(distance(gl,duppl)<pr && current_world==duppw && e2pbc(gl,duppl,duppw)){
                            float grad=1-distance(gl,duppl)/pr;
                            if(distance(gl,duppl2)<pr) FragColor=vec4(1,grad,1,1);
                            else FragColor=vec4(0,grad,1,1);
                        } else{
                            if(length(gl)<tpi&&iterations<10){
                                FragColor = vec4((gl.x+current_world)/3, (.5+T/2)*float((int(floor(gl.x))%2+int(floor(gl.y+current_world))%2+6)%2)/2+I/2+0.01,(gl.y)/3,1);
                            // FragColor=vec4(positions[66559],i,0,1);
                            } else FragColor=vec4(0.5,0.5,0.5,1);
                            FragColor*=vec4(ambientoc,ambientoc,ambientoc,1);
                        }
                    } else {
                        crl=sqrt(crl);
                        ccp=mat2(cp2.x,-cp2.y, cp2.y,cp2.x)*ccp/crl;
                        if(type==0) ccp/=crl;
                        //cd=1-cd/wallthickness/wallthickness;
                        //FragColor=vec4(cd*ccp[0],cd,cd*ccp[1],1);
                        FragColor=vec4(ccp[0],1,ccp[1],1);
                    }//*/
                    d=0;
                    //to do: add entity drawing
                }
            }//end E2
            if(iaunit > 0 && iterations < 9){//S2
                //d=2*atan(d);//test stereographic view
                bool wraparound = d > tpi;
                float dc = cos(d), ds = sin(d), rc = dot(sp1,sp2), rs = sqrt(1-rc*rc);
                mat3 rot;
                if(si < 0) rot = rotstart;
                else rot = s2matto(si);
                rot = lm*rot;
                if(si > -1){
                    rot = mat3(rc,0,rs, 0,1,0, -rs,0,rc)*rot;
                    rot = am*rot;
                }
                int i=int(positions[current_world]),cii=-1;
                vec3 cip1,cip2,cip;
                float cidr = s2disrank(dc,ds);
                if(wraparound) cidr = 5;//impossiblely high
                while(i<int(positions[current_world+1])){
                    vec3 p1=vec3(positions[i],positions[i+1],positions[i+2]), p2=vec3(positions[i+3],positions[i+4],positions[i+5]);
                    rc=dot(p1,p2);
                    rs=1-rc*rc;//actually rs^2, just reusing space
                    p1=rot*p1;
                    p2=rot*p2;
                    if(p1.y*p1.y<rs){
                        vec3 cipp, cipm;//closest intersection point plus and minus
                        float s = sqrt(rs - p1.y*p1.y), D = 1/(1-p1.y*p1.y);
                        cipp = vec3(p1.x*rc+p1.z*s,0,p1.z*rc-p1.x*s)*D;
                        cipm = vec3(p1.x*rc-p1.z*s,0,p1.z*rc+p1.x*s)*D;
                        float al = positions[i + 6];//angle limit, (cos(a)-1), decreases with angle
                        float ala = 1 + rs*al;
                        float disrank;
                        disrank=s2disrank(cipp.z,cipp.x);
                        if( (al==-2 || dot(cipp,p2)>ala ) && disrank<cidr && (si!=i || cipm.z>cipp.z)){//either not the source portal or the other point is closer
                            cii=i;
                            cip=cipp;
                            cidr=disrank;
                        }
                        disrank=s2disrank(cipm.z,cipm.x);
                        if( (al==-2 || dot(cipm,p2)>ala ) && disrank<cidr && (si!=i || cipp.z>cipm.z)){//either not the source portal or the other point is closer
                            cii=i;
                            cip=cipm;
                            cidr=disrank;
                        }
                        if(cii == i){
                            cip1=p1;
                            cip2=p2;
                        }
                    }
                    i += sizeOfpow;
                }
                if(wraparound && cii < 0){
                    cii = si;
                } else {
                    wraparound = false;
                }
                mat3 roti=transpose(rot);//inverse
                if(cii > -1){//hit something
                    int dat=int(positions[cii+7]);
                    if( (dat&4)>0 ) {
                        hitwall = true;
                        FragColor=vec4(0,0,0,1);
                        d=0;
                    } else {
                        float side=1;
                        if(((dat>>22)&1)==1) side=-1;
                        mirror=1;
                        if(((dat>>23)&1)==1) mirror=-1;
                        current_world=(dat>>3)&511;
                        int di = int(positions[current_world])+sizeOfpow*(dat>>12)&1023;
                        sp1 = vec3(positions[di],positions[di+1],positions[di+2]);
                        sp2 = vec3(positions[di+3],positions[di+4],positions[di+5]);
                        si = di;
                        //
                        if(wraparound){
                          d -= tpi;
                          am[0][1] *= mirror;
                          am[1][0] *= mirror;
                          am *= side;
                          am[2][2] *= side;
                          lm[0][1] *= mirror*side;
                          lm[1][0] *= mirror*side;
                          //...
                        } else {
                            d -= arctan(cip.x,cip.z);
                            float cosr2 = dot(cip1,cip2); cosr2 *= cosr2;
                            float c = (dot(cip2,cip) - cosr2)/(1 - cosr2);
                            float s = safe_sqrt(1 - c*c) * sign(dot(cip,cross(cip1,cip2)));
                            lm = mat3(c,-s*side*mirror,0, s*side*mirror,c,0, 0,0,1);
                            //lm = s2matto(cip1);//borrowing this memory space
                            //to do: find a faster way than s2matto
                            cip1 = mat3(cip.z,0,cip.x, 0,1,0, -cip.x,0,cip.z)*cip1;
                            //cip1 /= sqrt(cip1.x*cip1.x+cip1.y*cip1.y);
                            cip1 = normalize(cip1*vec3(1,1,0));
                            am = mat3(-cip1.x,mirror*-cip1.y,0, mirror*cip1.y,-cip1.x,0, 0,0,side)*side;
                            //cip2=lm*cip2;
                            //cip=lm*cip;
                            //cip=mat3(cip2.x,-cip2.y,0, cip2.y,cip2.x,0, 0,0,1)*cip;
                            //t=sqrt(cip.x*cip.x+cip.y*cip.y);
                            //cip/=t;
                            //lm=mat3(cip.x,-side*mirror*cip.y,0, side*mirror*cip.y,cip.x,0, 0,0,1);
                        }

                    }
                    /*

                    else if(distance(gl,duppl)<pr && current_world==duppw && e2pbc(gl,duppl,duppw)){
                        float grad=distance(gl,duppl)/pr;
                        if(distance(gl,duppl2)<pr) FragColor=vec4(1,grad,1,1);

                    */
                } else {
                    if(wraparound){
                        hitwall = true;
                        FragColor=vec4(0,0,0,1);
                    }
                    endOfRay = roti*vec3(ds,0,dc);
                    d=0;
                }
                //end of no portal hit
            }//end of s2
            if(iaunit < 0 && iterations < 9){//h2
                //d=2*atanh(d);
                float dc = exp(d), ds = 1/dc;
                dc = (dc + ds)*0.5; ds = dc - ds;
                //^cosh and sinh of d, done weird to remove redundant calculation
                float k = innerdot(sp1,sp2);
                mat3 rot;
                if(type > 0){
                    vec3 tempvec = sp1;
                    sp1 = sp2;
                    sp2 = tempvec;
                }
                rot = h2matto(sp1);
                sp2 = rot*sp2;
                t = length(sp2*vec3(1,1,0));
                sp2 /= t;
                rot = mat3(sp2.x,-sp2.y,0, sp2.y,sp2.x,0, 0,0,1)*rot;
                if(type == 2) rot = mat3(t,0,k, 0,1,0, k,0,t)*rot;
                rot = lm*rot;
                if(type == 2) rot = mat3(t,0,-k, 0,1,0, -k,0,t)*rot;
                if(type == 0 && si > 0) rot = mat3(-k,0,-t, 0,1,0, -t,0,-k)*rot;
                if(si > 0) rot = am*rot;
                vec3 cip1, cip2, cip;
                cip.x = ds;
                int cii = -1;
                int i = int(positions[current_world]);
                while(i < int(positions[current_world + 1])){
                    type = int(positions[i + 7])&3;
                    float limit = positions[i + 6];
                    vec3 p1 = vec3(positions[i],positions[i + 1],positions[i + 2]), p2 = vec3(positions[i + 3],positions[i + 4],positions[i + 5]);
                    k = innerdot(p1,p2);
                    p1 = rot*p1; p2 = rot*p2;
                    float o = sqrt(k*k - p1.z*p1.z + p1.x*p1.x);
                    float tR = (-k + o)/(p1.z - p1.x), itR = 1/tR;
                    float tL = (-k + o)/(p1.z + p1.x), itL = 1/tL;
                    vec3 cipR, cipL;
                    cipR.y = 0; cipL.y = 0;
                    cipR.x = (tR - itR)*0.5; cipR.z = cipR.x + itR;
                    cipL.x = -(tL - itL)*0.5; cipL.z = cipL.x + tL;
                    bool risvalid = tR > 0 && limit < innerdot(p2,cipR), lisvalid = tL > 0 && limit < innerdot(p2,cipL);
                    if(risvalid && 0 < cipR.x && cipR.x < cip.x && (i != si || (lisvalid && cipL.z < cipR.z) ) ){
                        cip = cipR;
                        cii = i;
                    }
                    if(lisvalid && 0 < cipL.x && cipL.x < cip.x && (i != si || (risvalid && cipL.z > cipR.z) ) ){
                        cip = cipL;
                        cii = i;
                    }
                    if(i == cii){
                        cip1 = p1;
                        cip2 = p2;
                    }
                    i += sizeOfpow;
                }//*/
                //...
                if(cii > 0){
                    int dat = int(positions[cii + 7]);
                    type = dat&3;
                    current_world=(dat>>3)&511;
                    int di = int(positions[current_world]) + sizeOfpow*(dat>>12)&1023;
                    sp1 = vec3(positions[di],positions[di + 1],positions[di + 2]);
                    sp2 = vec3(positions[di + 3],positions[di + 4],positions[di + 5]);
                    si = di;
                    float side = 1;
                    if(((dat>>22)&1)==1) side = -1;
                    mirror = 1;
                    if(((dat>>23)&1)==1) mirror = -1;
                    d -= acosh(cip.z);
                    k = innerdot(cip1,cip2);
                    float k2 = k*k;
                    float cL = -innerdot(cip,cip2), sL = 1;
                    if(type == 0){
                        cL = (cL - k2)/(1 - k2);
                        sL = 1 - cL*cL;
                    }
                    if(type == 1){//horocycle
                        sL = 2*cL - 2;
                    }
                    if(type == 2){//hypercycle
                        cL = (cL + k2)/(k2 + 1);
                        sL = cL*cL - 1;
                    }
                    sL = safe_sqrt(sL) * sign(innerdot(lcross(cip1,cip2),cip))*side*mirror;
                    if(type == 0){//circle
                        lm = mat3(cL,-sL,0, sL,cL,0, 0,0,1);
                    }
                    if(type == 1){//horocycle
                        lm = mat3(1 - sL*sL/2,-sL,-sL*sL/2, sL,1,sL,  sL*sL/2,sL,1 + sL*sL/2);
                    }
                    if(type == 2){//hypercycle
                        lm = mat3(1,0,0, 0,cL,sL, 0,sL,cL);
                    }
                    //calculate the sign and cos of the hit angle in a sense
                    cip1 = mat3(cip.z,0,-cip.x, 0,1,0, -cip.x,0,cip.z) * cip1;
                    if(type > 0) cip1 = -cip1;
                    cip1 = normalize(cip1*vec3(1,1,0));
                    am = mat3(-cip1.x,mirror*-cip1.y,0, mirror*cip1.y,-cip1.x,0, 0,0,side)*side;
                    //FragColor = vec4(Mod((cL),1),1,sL,1);
                    //FragColor = vec4(0*cL,1,sL,1);
                    //FragColor = vec4(-am[1][0],0,am[1][0],1);
                    //d = 0;hitwall = true;
                }
                else{
                //<temp>
                endOfRay = h2invert(rot)*vec3(ds,0,dc);
                d = safe_sqrt(endOfRay.z*endOfRay.z - 1);
                float x = endOfRay.x/d, y = endOfRay.y/d;
                float D = 1 - asinh(d)/2.5;
                D = sign(D);
                d = d/(endOfRay.z + 1);
                x = d*x, y = d*y;
                //^hyperboloid to poincare disk
                vec2 V = vec2(x*(1 - y) + x*(y + 1), 1 - y*y - x*x)/( (1 - y)*(1 - y) + x*x);
                //^disk to half plane
                float l = floorlog2(V.y);
                float R = Mod(l,2);
                float G = Mod(V.x/exponentialFloorb2(V.y),2);
                float B = Mod(V.x/exponentialFloorb2(V.y),1);
                if(G > 1){
                    G = 0;
                } else B = 0;
                FragColor = vec4(R,G,B,1/D)*max(D,0)+0.1;
                if(innerdot(endOfRay,pl) > -cosh(pr) && pw == current_world)
                    FragColor = vec4(0,1,1,1);
                d=0;
                /*float dfcc = 9;
                for(int n = 0; n < cellCount; n++){
                    dfcc = min(dfcc,-innerdot(endOfRay,cells[n]));
                }
                dfcc = 1/dfcc;
                //if(asinh(safe_sqrt(endOfRay.z*endOfRay.z - 1))<3)
                FragColor = vec4(dfcc,dfcc,dfcc,1);*/
                //FragColor *= dfcc*dfcc*dfcc*dfcc;
                //</temp>
                }
            }//end of h2
        }//end big while
        //
        if(!hitwall){
            int i = int(positions[current_world]);
            int nearstWallIndex = -1;
            float a = wallcornershad*wallcornershad;
            if(iaunit == 0){
                vec2 closestPoint, closestClosestPoint;
                float nwdis2 = wallcornershad*wallcornershad;
                while(i < positions[current_world + 1]){
                    int dat = int(positions[i + 7]);
                    if((dat&4) > 0){
                        vec2 p1 = vec2(positions[i],positions[i+1]),
                        p2 = vec2(positions[i+3],positions[i+4]) - p1,
                        endOfRayCopy = vec2(endOfRay.x,endOfRay.y) - p1;
                        float rl2 = dot(p2,p2);
                        int type = dat&3;
                        if(type == 0){
                            closestPoint = endOfRayCopy*sqrt(rl2/dot(endOfRayCopy,endOfRayCopy));;
                            float al = positions[i + 6], ala = al*rl2,
                            dis = distance(closestPoint, p2);
                            dis *= dis;
                            if(al < 4 && dis > ala){
                                /*vec2 ep1 = vec2(positions[i + 8],positions[i + 9]) - p1;
                                float d1 = distance2(ep1, endOfRayCopy);
                                vec2 ep2 = vec2(positions[i + 11],positions[i + 12]) - p1;
                                float d2 = distance2(ep2, endOfRayCopy);
                                if(d1 < d2) closestPoint = ep1;
                                else closestPoint = ep2;*/
                                float c = 1 - al*0.5, s = sqrt(1 - c*c);
                                if(endOfRayCopy.y*p2.x - endOfRayCopy.x*p2.y < 0) s = -s;
                                closestPoint = mat2(c,s, -s,c)*p2;
                            }
                        } else {
                            closestPoint = p2*min(max(0,dot(endOfRayCopy,p2)),rl2)/rl2;
                        }
                        float dis2 = distance(closestPoint, endOfRayCopy);
                        dis2 *= dis2;
                        if(dis2 < wallcornershad*wallcornershad) a -= dis2;
                        if(dis2 < 0.99*nwdis2){
                            nwdis2 = dis2;
                            nearstWallIndex = i;
                            closestClosestPoint = closestPoint;
                        }
                    }
                    i += sizeOfpow;
                }
                if(distance(endOfRay,pl) < pr && current_world == pw && !e2pbc(endOfRay,pl,pw)){
                    float grad=1-distance(endOfRay,pl)/pr;
                    if(distance(endOfRay,pl2) < pr) FragColor=vec4(1,grad,1,1);
                    else FragColor=vec4(0,grad,1,1);
                } else if(distance(endOfRay,duppl) < pr && current_world == duppw && e2pbc(endOfRay,duppl,duppw)){
                    float grad=1-distance(endOfRay,duppl)/pr;
                    if(distance(endOfRay,duppl2) < pr) FragColor=vec4(1,grad,1,1);
                    else FragColor=vec4(0,grad,1,1);
                } else {
                    vec2 z = vec2(endOfRay.x,endOfRay.y);
                    int N = 30, n = N;
                    while(n > 0 && dot(z,z) < 4){
                        t = z.x*z.x - z.y*z.y + endOfRay.x;
                        z.y = 2*z.x*z.y + endOfRay.y;
                        z.x = t;
                        n--;
                    }
                    I=float(n)/float(N);
                    float test = positions[66559];
                    if(dot(endOfRay,endOfRay) < tpi*tpi){
                        FragColor = vec4((endOfRay.x+current_world)/3, (.5+test/2)*float((int(floor(endOfRay.x))%2+int(floor(endOfRay.y+current_world))%2+6)%2)/2+I/2+0.01,(endOfRay.y)/3,1);
                    }else FragColor=vec4(0.5,0.5,0.5,1);
                    if(nearstWallIndex > 0){
                        i = nearstWallIndex;
                        vec2 p1 = vec2(positions[i],positions[i+1]),
                        p2 = vec2(positions[i+3],positions[i+4]) - p1,
                        endOfRayCopy = vec2(endOfRay.x,endOfRay.y) - p1;
                        if(nwdis2 < wallthickness*wallthickness){
                            float relLoc;
                            int type = int(positions[i + 7])&3;
                            relLoc = distance(closestClosestPoint,p2);
                            relLoc *= relLoc;
                            if(type == 0){
                                //if performance allows:
                                //float rr = dot(p2,p2);
                                //relLoc = acos(1 - relLoc/rr*0.5) * sqrt(rr);
                                if(dot(vec2(-p2.y,p2.x),closestClosestPoint) < 0) relLoc = -relLoc;
                                //...
                            } else {
                                relLoc = sqrt(relLoc);
                            }
                            FragColor = vec4(relLoc,-relLoc,1,1);
                        } else {
                            a = sqrt(nwdis2)/wallcornershad;
                            FragColor *= vec4(a,a,a,1);
                            //...darken floor here
                        }
                        //...
                    }
                }
                //...
            }
            if(iaunit > 0){
                    vec3 closestPoint, closestClosestPoint;
                    float codtccp = coswallcornershad;//cos of distance to closest closest point, I am not typing that all out outside of this comment
                    while(i < positions[current_world + 1]){
                        if( (int(positions[i + 7])&4) > 0){
                            vec3 p1 = vec3(positions[i],positions[i + 1],positions[i + 2]),
                            p2 = vec3(positions[i + 3],positions[i + 4],positions[i + 5]);
                            float rc = dot(p1,p2), rss = 1 - rc*rc;//cos of r and sin(r)^2
                            float Retro = dot(endOfRay,p1),
                            A = safe_sqrt( rss/(1 - Retro*Retro ) ),
                            B = rc - A * Retro;
                            closestPoint = A*endOfRay + B*p1;
                            float al = positions[i+6], ala = 1 + al*rss;
                            if(al>-2 && dot(closestPoint,p2) < ala){
                                vec3 ep1 = vec3(positions[i + 8],positions[i + 9],positions[i + 10]);
                                vec3 ep2 = vec3(positions[i + 11],positions[i + 12],positions[i + 13]);
                                if(dot(ep1,endOfRay) > dot(ep2,endOfRay)) closestPoint = ep1;
                                else closestPoint = ep2;
                            }
                            float cosDisToCP = dot(endOfRay, closestPoint);
                            if(cosDisToCP > codtccp*1.00001 && Retro < pr){//no "z" fighting and bounds checking
                                nearstWallIndex = i;
                                closestClosestPoint = closestPoint;
                                codtccp = cosDisToCP;
                            }
                        }
                        i += sizeOfpow;
                    }
                    if(dot(endOfRay,pl)>prc && current_world==pw && !s2pbc(endOfRay,pl,current_world)){
                        float grad = 1 - smallArccos(dot(endOfRay,pl))/pr;
                        if(dot(endOfRay,pl2) > prc) FragColor = vec4(1,grad,1,1);
                        else FragColor=vec4(0,grad,1,1);
                    } else if(dot(endOfRay,duppl)>prc && current_world==duppw && s2pbc(endOfRay,duppl,current_world)){
                        float grad = 1 - smallArccos(dot(endOfRay,duppl))/pr;
                        if(dot(endOfRay,duppl2)>prc) FragColor = vec4(1,grad,1,1);
                        else FragColor = vec4(0,grad,1,1);
                    } else {
                        float k=sign(endOfRay.x*endOfRay.y*endOfRay.z);
                        vec3 color = (sign(endOfRay) + 1 - endOfRay)/2;
                        FragColor = vec4(color.x,color.y,color.z,1);
                        if(nearstWallIndex>0) {
                            if(codtccp > coswallthickness){
                            vec3 p1 = vec3(positions[nearstWallIndex],positions[nearstWallIndex + 1],positions[nearstWallIndex + 2]),
                            p2 = vec3(positions[nearstWallIndex + 3],positions[nearstWallIndex + 4],positions[nearstWallIndex + 5]);
                            float relLoc = dot(endOfRay,p2) * sign(dot(cross(p1,p2),endOfRay));
                            FragColor = vec4(relLoc,1,-relLoc,1);
                        } else {
                                k = (safe_sqrt( 2 - 2*codtccp))/(wallcornershad);
                                FragColor *= vec4(k,k,k,1);
                        }
                    }
                }
            }
        }
        if(iterations>8) FragColor = bgcolor;
        if(iterations<1) FragColor=vec4(1,1,1,1);//banadaid
    }
    else {//outside the view area
        FragColor = bgcolor;//8
    }
    //FragColor=round(FragColor);
}
