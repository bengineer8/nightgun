#version 460 core

#define portsAndWallsSize 1048576
#define maxits 9


#define datWorldOff 3
#define datWorldMask 2047
#define datIndexOff 14
#define datIndexMask 255

bool debug = false;
const float tpi=6.283185307;//two pi
//const int shaderDataLength=8*128+8*4096+4*8192;
const int sizeOfpow = 14;
uniform float pr = 1.0/9.0;
const float wallthickness=1/36.0;
const float wallcornershad=1.0/9;
const float coswallthickness=cos(wallthickness);
const float coswallcornershad=cos(wallcornershad);
const float coshwallthickness=cosh(wallthickness);
const float coshwallcornershad=cosh(wallcornershad);
const float inf = 1.0/0;
const mat3 negninty = mat3(0,-1,0, 1,0,0, 0,0,1);
const mat3 ninty = mat3(0,1,0, -1,0,0, 0,0,1);
//const vec4 bgcolor=vec4(0,.125,0.498,1);
const vec4 bgcolor=vec4(0,0,0,1);
out vec4 FragColor;
uniform int res=800;
uniform float T=1.0,zoom=1.0;
//uniform vec2 playerFacing = vec2(1,0);
uniform vec3 pl,pl2,camRef,playerduploc,playerduploc2;
uniform int pw;//player world
uniform int duppw;
uniform float mirrorPlayer = 1;//-1 for mirror player
uniform float mirrorDup = 1;
uniform int iopip;
uniform int diopip;

layout (std430, binding=3) buffer ssbo1
{
  float portsAndWalls[portsAndWallsSize];
};
uniform float iaunits[128];

void disSqrstoxy(float ds1, float ds2, float r, float s, out float x, out float y){
    y = (ds1 - ds2 + r*r)/r*0.5;
    x = s*sqrt(ds1 - y*y);
}

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

float lidot(vec3 a, vec3 b){
    return dot(a,b*vec3(-1,-1,1));
}

vec3 icross(vec3 a, vec3 b){
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

float safe_acos(float x){
    return acos(min(1,max(-1,x)));
}

mat3 s2matto(vec3 p){//returns a matrix that rotates a unit sphere so p is at (0,0,1)
    float r=sqrt(p.x*p.x+p.y*p.y), c=p.x/r, s=p.y/r;
    if(r==0){
        c=1;s=0;
    }
    return mat3(p.z*c,-s,p.x, p.z*s,c,p.y, -r,0,p.z);
}

mat3 s2matto(int i){
    float x = portsAndWalls[i],
    y = portsAndWalls[i + 1],
    z = portsAndWalls[i + 2],
    ns = portsAndWalls[i + 8],
    c = portsAndWalls[i + 9],
    nr = portsAndWalls[i + 10],
    c2 = portsAndWalls[i + 11],
    s2 = portsAndWalls[i + 12];
    return mat3(c2,-s2,0, s2,c2,0, 0,0,1)*mat3(z*c,ns,x, -z*ns,c,y, nr,0,z);
    //return mat3(1,0,0, 0,1,0, 0,0,1);
}

bool s2PortalBetween(vec3 A, vec3 B, int i){
    if(i > 0){
        vec3 p1 = {portsAndWalls[i], portsAndWalls[i + 1], portsAndWalls[i + 2]}, p2 = {portsAndWalls[i + 3], portsAndWalls[i + 4], portsAndWalls[i + 5]};
        float k = dot(p1,p2);
        if((dot(A,p1) - k)*(dot(B,p1) - k) < 0) return true;
    }
    return false;
}


bool e2PortalBetween(vec3 A, vec3 B, int i){
    if(i > 0){
        vec3 p1 = {portsAndWalls[i], portsAndWalls[i + 1], portsAndWalls[i + 2]}, p2 = {portsAndWalls[i + 3], portsAndWalls[i + 4], portsAndWalls[i + 5]};
        int type = int(portsAndWalls[i + 7])&1;
        float k,a,b;
        if(type == 0){
            k = distance(p1,p2);
            a = distance(A,p1);
            b = distance(B,p1);
        } else{
            p2 = mat3(0,1,0, -1,0,0, 0,0,1)*(p2 - p1);
            k = dot(p1,p2);
            a = dot(A,p2);
            b = dot(B,p2);
        }
        if((a - k)*(b - k) < 0) return true;

    }
    return false;
}

float e2side(vec3 a, vec3 b, vec3 c){
    b = mat3(0,-1,0, 1,0,0, 0,0,1)*(b - a);
    c -= a;
    return sign(dot(b,c));
}

mat3 h2matto(vec3 p){
    float r = sqrt(p.x*p.x + p.y*p.y), c = p.x/r, s = p.y/r;
    if(r == 0){
        c = 1;
        s = 0;
    }
    return mat3(p.z*c,-s,-p.x, p.z*s,c,-p.y, -r,0,p.z);
}

bool h2PortalBetween(vec3 A, vec3 B, int i){
    if(i > 0){
        vec3 p1 = {portsAndWalls[i], portsAndWalls[i + 1], portsAndWalls[i + 2]}, p2 = {portsAndWalls[i + 3], portsAndWalls[i + 4], portsAndWalls[i + 5]};
        float k = lidot(p1,p2);
        if((lidot(A,p1) - k)*(lidot(B,p1) - k) < 0) return true;
    }
    return false;
}

void main(){
    FragColor = vec4(0,0,0,1);
    vec4 dcolor = vec4(0,0,0,1);
    int previousPortal = -1;
    int previousPreviousPortal = -1;
    bool ignoreNextSelf = false;
    int current_world=pw;
    vec2 screenPos=vec2(2.0*gl_FragCoord.x/float(res)-1.0,2.0*gl_FragCoord.y/float(res)-1.0);
    float r,t;
    float mirror=mirrorPlayer;
    if(length(screenPos)<1.0){
        screenPos.y *= mirrorPlayer;
        int type = -1, sourceIndex = -1;
        float d = length(screenPos);
        screenPos = screenPos/d;
        d = d*tpi*zoom;
        vec3 sourcept1 = pl, sourcept2 = camRef;
        mat3 entranceAng = mat3(screenPos[0],-screenPos[1],0,screenPos[1],screenPos[0],0,0,0,1);
        float entranceLoc;
        vec3 endOfRay;
        float iaunit;
        bool hitwall = false;
        int iterations = 0;
        while(d>0 && iterations < maxits){
            iaunit = iaunits[current_world];
            //...
            if(iaunit < 0){//BUG arc portals no fully worky
                iaunit *= -1;
                float k = lidot(sourcept1,sourcept2);
                float ds = sinh(d*iaunit), dc = sqrt(1 + ds*ds);
                if(type > 0){
                    vec3 tv = sourcept1;
                    sourcept1 = sourcept2;
                    sourcept2 = tv;
                }
                mat3 transform = h2matto(sourcept1);
                sourcept2 = transform*sourcept2;
                //float k = sourcept2.z;
                r = length(sourcept2 * vec3(1,1,0));
                //r = sqrt(k*k - 1 + max(0,type));
                sourcept2 /= r;
                if(type > 0){
                    sourcept2 *= -1;
                }
                transform = mat3(sourcept2.x,-sourcept2.y,0, sourcept2.y,sourcept2.x,0, 0,0,1)*transform;
                //r = sqrt(k*k - 1 + max(0,type));
                entranceLoc *= -iaunit/r;
                if(type == 0){
                    float c = cos(entranceLoc), s = sin(entranceLoc);
                    transform = mat3(c,s,0, -s,c,0, 0,0,1)*transform;
                    transform = mat3(k,0,-r, 0,1,0, -r,0,k)*transform;
                } else if(type == 1){
                    entranceLoc *= r;
                    float A = entranceLoc*entranceLoc*0.5;
                    transform = mat3(1 - A,entranceLoc,A, -entranceLoc,1,entranceLoc, -A,entranceLoc,1 + A)*transform;
                } else if(type == 2){
                    transform = mat3(r,0,k, 0,1,0, k,0,r)*transform;
                    float s = sinh(entranceLoc), c = sqrt(s*s + 1);
                    transform = mat3(1,0,0, 0,c,s, 0,s,c)*transform;
                    transform = mat3(r,0,-k, 0,1,0, -k,0,r)*transform;
                }
                transform = entranceAng*transform;
                vec3 cip1, cip2, cip;
                cip.x = ds;
                int cii = -1;
                int i = int(portsAndWalls[current_world]);
                float cik;
                while(i < int(portsAndWalls[current_world + 1])){
                    float limit = portsAndWalls[i + 6];
                    vec3 p1 = vec3(portsAndWalls[i],portsAndWalls[i + 1],portsAndWalls[i + 2]), p2 = vec3(portsAndWalls[i + 3],portsAndWalls[i + 4],portsAndWalls[i + 5]);
                    k = lidot(p1,p2);
                    p1 = transform*p1; p2 = transform*p2;
                    float o = sqrt(k*k - p1.z*p1.z + p1.x*p1.x);
                    if(p1.z<0) o = -o;
                    o += k;
                    float tR = o/(p1.z + p1.x);
                    float tL = (p1.z - p1.x)/o;
                    vec3 cipR, cipL;
                    cipR.y = 0; cipL.y = 0;
                    cipR.x = (1/tR - tR)*0.5; cipR.z = cipR.x + tR;
                    cipL.x = (1/tL - tL)*0.5; cipL.z = cipL.x + tL;
                    bool risvalid = tR >= 0 && limit > lidot(p2,cipR), lisvalid = tL >= 0 && limit > lidot(p2,cipL);
                    if(!(sourceIndex == i && ignoreNextSelf)){
                        if(risvalid && 0 < cipR.x && cipR.x < cip.x && (i != sourceIndex || (lisvalid && cipL.z < cipR.z) ) ){
                            cip = cipR;
                            cii = i;
                            }
                        if(lisvalid && 0 < cipL.x && cipL.x < cip.x && (i != sourceIndex || (risvalid && cipL.z > cipR.z) ) ){
                            cip = cipL;
                            cii = i;
                        }
                        if(!(o == o && o != 0) && i == sourceIndex && previousPreviousPortal != i){//this can only be caused by rounding errors as it should at least detect the point it came from
                            //if it hit the portal 2 iterations ago^, that means it is trying to exit and thus an intersection should't be forced
                            cii = i;
                            cip = vec3(0,0,1);
                            ignoreNextSelf = true;
                        }//^hack fix
                        if(i == cii){
                            cip1 = p1;
                            cip2 = p2;
                            cik = k;
                        }
                    } else ignoreNextSelf = false;
                    i += sizeOfpow;
                }//*/
                if(cii > 0){
                    int dat = int(portsAndWalls[cii + 7]);
                    type = dat&7;
                    if(type > 3){
                        d = 0;
                        hitwall = true;
                    } else {
                        previousPreviousPortal = previousPortal;
                        previousPortal = cii;
                        current_world = (dat>>datWorldOff)&datWorldMask;
                        int di = int(portsAndWalls[current_world]) + sizeOfpow*((dat>>datIndexOff)&datIndexMask);
                        sourcept1 = vec3(portsAndWalls[di],portsAndWalls[di + 1],portsAndWalls[di + 2]);
                        sourcept2 = vec3(portsAndWalls[di + 3],portsAndWalls[di + 4],portsAndWalls[di + 5]);
                        sourceIndex = di;
                        float side = 1;
                        if(((dat>>22)&1)==1) side = -1;
                        mirror = 1;
                        if(((dat>>23)&1)==1) mirror = -1;
                        //
                        float aunit = 1/iaunit;
                        cip2.z = sqrt(cip2.x*cip2.x + cip2.y*cip2.y + 1);//correcting rounding
                        entranceLoc = lidot(cip,cip2);
                        //k = lidot(cip1,cip2);
                        k = cik;
                        float updown = sign(lidot(cip,icross(cip2,cip1)));
                        cip1 = mat3(cip.z,0,-cip.x, 0,1,0, -cip.x,0,cip.z)*cip1;
                        r = k*k - 1 + type;
                        if(type == 0){
                            entranceLoc = safe_acos((k*k - entranceLoc)/r);//I have a backup method if this continues to b a problem after chunking
                        } else if(type == 1){
                            entranceLoc = safe_sqrt(entranceLoc*2 - 2)/k;
                        } else {
                            entranceLoc = acosh(max((k*k + entranceLoc)/r,1));
                        }
                        r = safe_sqrt(r)*side;
                        entranceLoc = r*entranceLoc*updown*mirror*aunit;
                        cip1 = -normalize(cip1*vec3(1,mirror,0))*side;//NOTE: On paper, I would only need to divide by the r above, but, rounding says "no"
                        entranceAng = mat3(cip1.x,cip1.y,0, -cip1.y,cip1.x,0, 0,0,1);
                        type = int(portsAndWalls[di + 7])&3;
                        d -= aunit*acosh(cip.z);
                        //FragColor = vec4(entranceLoc,0,-entranceLoc,1); hitwall = true; d = 0;
                    }
                } else {
                    mat3 itransform = h2invert(transform);
                    endOfRay = itransform*vec3(ds,0,dc);
                    d = 0;
                }
                //...
                iaunit *= -1;
            } if(iaunit == 0){
                vec3 off = -sourcept1;
                sourcept2 += off;
                r = length(sourcept2);
                sourcept2 /= r;
                mat3 rot = mat3(sourcept2.x,-sourcept2.y,0, sourcept2.y,sourcept2.x,0, 0,0,1);
                off = rot*off;
                if(type == 0){
                    float a = entranceLoc/r;
                    float c = cos(a), s = -sin(a);
                    mat3 tm = mat3(c,s,0, -s,c,0, 0,0,1);
                    off = tm*off;
                    rot = tm*rot;
                    off.x -= r;
                }
                if(type == 1){
                    off.x -= entranceLoc + r/2;
                    off = ninty*off;
                    rot = ninty*rot;
                }
                off = entranceAng*off;
                rot = entranceAng*rot;
                int i = int(portsAndWalls[current_world]),cii = -1;
                vec3 cip1, cip2;
                float cid = d;
                while(i<int(portsAndWalls[current_world+1]) && iterations<maxits){
                    vec3 p1=vec3(portsAndWalls[i],portsAndWalls[i+1],0), p2=vec3(portsAndWalls[i+3],portsAndWalls[i+4],0);
                    t = distance(p1,p2);
                    type = int(portsAndWalls[i + 7])&3;
                    float limit = portsAndWalls[i + 6];
                    p1 = rot*p1 + off; p2 = rot*p2 + off;
                    if(!(sourceIndex == i && ignoreNextSelf)){
                        if(type == 0){
                            float o = sqrt(t*t - p1.y*p1.y);
                            if(abs(p1.y)<= t && (p2.y*p2.y <= limit)){
                                float xt = p1.x - o, xt1 = xt;
                                if(((xt - p2.x)*(xt - p2.x) + p2.y*p2.y < limit) && xt >= 0 && xt<cid && (i != sourceIndex)) {//1st collision is not allowed if portal is the source
                                    cii=i; cid = xt;
                                }
                                xt = p1.x + o;
                                if(((xt - p2.x)*(xt - p2.x) + p2.y*p2.y < limit) && xt >= 0 && xt<cid && (i != sourceIndex||abs(xt1)<xt)) {//2nd clollision is only allowed with oneself if the 1st collsion would be the starting point
                                    cii=i; cid = xt;
                                }
                            }
                            if(!(o == o && o != 0) && i == sourceIndex && previousPreviousPortal != i){//Forcing things when a rouding error is detected
                                cid = 0;
                                cii = i;
                                ignoreNextSelf = true;
                            }
                        }
                        if(type == 1 && i != sourceIndex && p1.y*p2.y<=0 && (p1.x>0 || p2.x>0) ){
                            float xt=(p1.x-p2.x)*p2.y/(p2.y-p1.y)+p2.x;
                            if(xt>0 && xt<cid){
                                cii=i; cid=xt;
                            }
                        }
                    } else ignoreNextSelf = false;
                    if(cii == i){
                        cip1 = p1; cip2 = p2;
                    }
                    i += sizeOfpow;
                }
                if(cii > 0){
                    int dat=int(portsAndWalls[cii+7]);
                    type = dat&7;
                    if(type > 3){
                        hitwall = true;
                        d = 0;
                    } else {
                        cip1.x -= cid;
                        cip2.x -= cid;
                        r = length(cip1);
                        cip2 -= cip1;
                        cip1 *= -1;
                        entranceAng[0] = cip1/r;
                        float R = length(cip2);
                        float R2 = R*R;
                        if(type == 0){
                            entranceLoc = R*sign(cip1.y*cip2.x - cip1.x*cip2.y)*safe_acos(1 - dot(cip1-cip2,cip1-cip2)/R2*.5);
                        } else {
                            entranceAng[0] = negninty*entranceAng[0];
                            entranceLoc = r - R/2;
                            //...
                        }
                        float side = 1;
                        if(((dat>>22)&1)==1) side=-1;
                        mirror = 1;
                        if(((dat>>23)&1)==1) mirror=-1;
                        previousPreviousPortal = previousPortal;
                        previousPortal = cii;
                        type=dat&3;
                        current_world = (dat>>datWorldOff)&datWorldMask;
                        int di = int(portsAndWalls[current_world]) + sizeOfpow*((dat>>datIndexOff)&datIndexMask);
                        sourcept1 = vec3(portsAndWalls[di],portsAndWalls[di+1],portsAndWalls[di+2]);
                        sourcept2 = vec3(portsAndWalls[di+3],portsAndWalls[di+4],portsAndWalls[di+5]);
                        sourceIndex = di;
                        type = int(portsAndWalls[di + 7])&7;
                        //
                        entranceAng[0] *= side;
                        entranceAng[0].y *= mirror;
                        entranceLoc *= side*mirror;
                        entranceAng = mat3(entranceAng[0].x,entranceAng[0].y,0, -entranceAng[0].y,entranceAng[0].x,0, 0,0,1);
                        //d=0; hitwall = true; FragColor = vec4(entranceLoc,0,-entranceLoc,1);
                        //...
                        d -= cid;
                    }
                    //...
                } else {
                    mat3 roti = transpose(rot);
                    endOfRay = roti*(vec3(d,0,0) - off);
                    d = 0;
                }
                //...
            } if(iaunit > 0){
                float t = d*iaunit;
                bool wraparound = t > tpi;
                float dc = cos(t), ds = sin(t), rc = dot(sourcept1,sourcept2), rs = sqrt(1 - rc*rc);
                mat3 rot = s2matto(sourcept1);
                sourcept2 = rot*sourcept2/rs;
                rot = mat3(sourcept2.x,-sourcept2.y,0, sourcept2.y,sourcept2.x,0, 0,0,1)*rot;
                if(sourceIndex > -1){
                    float a = entranceLoc*iaunit/rs;
                    float c = cos(a), s = sin(a);
                    rot = mat3(rc,0,rs, 0,1,0, -rs,0,rc)*mat3(c,-s,0, s,c,0, 0,0,1)*rot;
                }
                rot = entranceAng*rot;
                int i = int(portsAndWalls[current_world]), cii = -1;
                vec3 cip1,cip2,cip;
                float cidr = s2disrank(dc,ds);
                if(wraparound) cidr = 5;//impossiblely high
                while(i<int(portsAndWalls[current_world+1])){
                    vec3 p1=vec3(portsAndWalls[i],portsAndWalls[i+1],portsAndWalls[i+2]), p2=vec3(portsAndWalls[i+3],portsAndWalls[i+4],portsAndWalls[i+5]);
                    rc=dot(p1,p2);
                    rs=1-rc*rc;//actually rs^2, just reusing space
                    p1=rot*p1;
                    p2=rot*p2;
                    if(!(sourceIndex == i && ignoreNextSelf)){
                        if(p1.y*p1.y <= rs){
                            vec3 cipp, cipm;//closest intersection point plus and minus
                            float s = sqrt(rs - p1.y*p1.y), D = 1/(1-p1.y*p1.y);
                            cipp = vec3(p1.x*rc+p1.z*s,0,p1.z*rc-p1.x*s)*D;
                            cipm = vec3(p1.x*rc-p1.z*s,0,p1.z*rc+p1.x*s)*D;
                            float limit = portsAndWalls[i + 6];
                            float disrank;
                            disrank=s2disrank(cipp.z,cipp.x);
                            if( dot(cipp,p2) > limit && disrank<cidr && (sourceIndex != i || cipm.z > cipp.z)){//either not the source portal or the other point is closer
                                cii = i;
                                cip = cipp;
                                cidr = disrank;
                            }
                            disrank=s2disrank(cipm.z,cipm.x);
                            if( dot(cipm,p2) > limit && disrank<cidr && (sourceIndex != i || cipp.z > cipm.z)){//either not the source portal or the other point is closer
                                cii = i;
                                cip = cipm;
                                cidr = disrank;
                            }
                            if(!(s == s && cipm.z != cipp.z) && i == sourceIndex && previousPreviousPortal != i){//Forcing things when a rouding error or knife edge is detected
                                cii = i;
                                cip = vec3(0,0,1);
                                cidr = 0;
                                ignoreNextSelf = true;
                            }
                            if(cii == i){
                                cip1=p1;
                                cip2=p2;
                            }
                        }
                    } else ignoreNextSelf = false;
                    i += sizeOfpow;
                }
                if(wraparound && cii < 0){
                    cii = sourceIndex;
                } else {
                    wraparound = false;
                }
                mat3 irot = transpose(rot);
                if(cii > 0){
                    int dat=int(portsAndWalls[cii+7]);
                    type = dat&7;
                    if(type > 3){
                        d = 0;
                        hitwall = true;
                    } else {
                        previousPreviousPortal = previousPortal;
                        previousPortal = cii;
                        int di;
                        float side = 1;
                        if(((dat>>22)&1)==1) side=-1;
                        mirror = 1;
                        if(((dat>>23)&1)==1) mirror=-1;
                        if(wraparound){
                            d -= tpi/iaunit;
                            di = cii;
                        } else {
                            current_world = (dat>>datWorldOff)&datWorldMask;
                            di = int(portsAndWalls[current_world]) + sizeOfpow*((dat>>datIndexOff)&datIndexMask);
                            sourceIndex = di;
                            type = int(portsAndWalls[di + 7])&7;
                            rc = dot(cip1,cip2); rs = 1 - rc*rc;
                            entranceLoc = dot(cip,cip2);
                            entranceLoc = safe_acos((entranceLoc - rc*rc)/rs)/iaunit*sign(dot(cip,cross(cip1,cip2)));
                            rs = sqrt(rs);
                            entranceLoc *= rs;
                            cip1 = mat3(cip.z,0,cip.x, 0,1,0, -cip.x,0,cip.z)*cip1/-rs;
                            cip1 = normalize(cip1*vec3(1,1,0));
                            entranceAng = mat3(cip1.x,cip1.y,0, -cip1.y,cip1.x,0, 0,0,1);
                            d -= arctan(cip.x,cip.z)/iaunit;
                        }

                        sourcept1 = vec3(portsAndWalls[di],portsAndWalls[di+1],portsAndWalls[di+2]);
                        sourcept2 = vec3(portsAndWalls[di+3],portsAndWalls[di+4],portsAndWalls[di+5]);
                        entranceAng *= side;
                        entranceAng[2][2] *= side;
                        entranceAng[0][1] *= mirror;
                        entranceAng[1][0] *= mirror;
                        entranceLoc *= mirror*side;
                    }
                    //d = 0;hitwall = true;
                } else {
                    endOfRay = irot * vec3(ds,0,dc);
                    d = 0;
                }
                //...
            }
            iterations++;
            //...
        }
        if(!hitwall && iterations <= maxits){
            int i = int(portsAndWalls[current_world]);
            int nearstWallIndex = -1;
            float a = wallcornershad*wallcornershad;
            //...
            if(iaunit < 0){
                vec2 V = vec2(endOfRay.x,endOfRay.y)/(endOfRay.z + 1);
                //^hyperboloid to poincare disk
                V = vec2(2*V.x, 1 - V.y*V.y - V.x*V.x)/( (1 - V.y)*(1 - V.y) + V.x*V.x);
                //^disk to half plane*/


                //V = vec2(endOfRay.x,1)*(endOfRay.z + 1)/(endOfRay.z*(endOfRay.z + 1 - endOfRay.y) - endOfRay.y);

                /*vec2 V;
                V.y = 1/(endOfRay.z - endOfRay.y);
                V.x = endOfRay.x*V.y;
                //hyperboloid to half plane^*/
                //^NOTE:errors creep in too fast to use

                float l = floorlog2(V.y);
                float R = Mod(l,2);
                float G = Mod(V.x/exponentialFloorb2(V.y),2);
                float B = Mod(V.x/exponentialFloorb2(V.y),1);
                if(G > 1){
                    G = 0;
                } else B = 0;
                FragColor = vec4(R,G,B,1);
                if(endOfRay.z > cosh(3)) FragColor *= vec4(0.5,0.5,0.5,1);//about where rounding errors are noticable
                //^debug tiling
                //
                float scale = 2;
                r = sqrt(2*abs(endOfRay.y));
                float A = sign(endOfRay.x)*asinh(abs(endOfRay.x)/sqrt(endOfRay.y*endOfRay.y + 1));
                vec2 P = scale*r*vec2(-cos(A)*sign(endOfRay.y),sin(A)) + vec2(1 + scale*sqrt(2),0)*sign(endOfRay.y);
                vec2 z = P;
                    int N = 30, n = N;
                    while(n > 0 && dot(z,z) < 4){
                        t = z.x*z.x - z.y*z.y + P.x;
                        z.y = 2*z.x*z.y + P.y;
                        z.x = t;
                        n--;
                    }
                float I = float(n)/float(N);
                //FragColor = vec4((mod(P.x-tpi,2*tpi)-tpi)/6,0.75*Mod(floor(P.x)+floor(P.y),2)/2+I/2+0.01,(mod(P.y-tpi,2*tpi)-tpi)/6,1);
                //FragColor = vec4((mod(P.x-tpi,2*tpi)-tpi)/6, (.75)*float((int(floor(P.x))%2+int(floor(P.y))%2+6)%2)/2+1/2,(mod(P.y-tpi,2*tpi)-tpi)/6,1);
                if(current_world == 32) FragColor = vec4((mod(P.x-tpi,2*tpi)-tpi)/6, 0.75*Mod(floor(P.x) + floor(P.y),2)/2+I/2+0.01,(mod(P.y-tpi,2*tpi)-tpi)/6,1);
                //...
                vec3 closestPoint;
                float coshdis = cosh(wallcornershad*-iaunit);
                while(i < portsAndWalls[current_world + 1]){
                    type = int(portsAndWalls[i + 7])&7;
                    if(type > 3){
                        vec3 p1 = vec3(portsAndWalls[i],portsAndWalls[i + 1],portsAndWalls[i + 2]),
                        p2 = vec3(portsAndWalls[i + 3],portsAndWalls[i + 4],portsAndWalls[i + 5]);
                        float
                        limit = portsAndWalls[i + 6],
                        k = lidot(p1,p2),
                        id1 = lidot(endOfRay,p1),
                        id2 = lidot(p1,p1),
                        A = id1*id1 - id2,
                        B = sqrt( (k*k - id2)/A ),
                        C;
                        if(type == 5) C = (id1*id1 - k*k)/(2*A*k);
                        else C = (k - B*id1)/id2;
                        vec3 cp = B*endOfRay + C*p1;
                        if(lidot(cp,p2) > limit){
                            int off;
                            if(lidot(icross(p1,p2),cp) > 0) off = 8;
                            else off = 11;
                            cp = vec3(portsAndWalls[i + off], portsAndWalls[i + off + 1], portsAndWalls[i + off + 2]);
                        }
                        float cd = lidot(endOfRay,cp);
                        if(cd < 0.99999*coshdis){
                            coshdis = cd;
                            nearstWallIndex = i;
                            closestPoint = cp;
                        }
                        //...
                    }
                    //...
                    i += sizeOfpow;
                }
                if(nearstWallIndex > 0){
                    i = nearstWallIndex;
                    vec3 p1 = vec3(portsAndWalls[i],portsAndWalls[i + 1],portsAndWalls[i + 2]),
                    p2 = vec3(portsAndWalls[i + 3],portsAndWalls[i + 4],portsAndWalls[i + 5]);
                    float k = lidot(p1,p2);
                    type = int(portsAndWalls[i + 7])&3;
                    float coshwallthickness = cosh(wallthickness*-iaunit);
                    if(coshdis < coshwallthickness){
                        float k2 = k*k;
                        float cL = lidot(closestPoint,p2), sL = 1;
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
                        //jank V
                        sL = safe_sqrt(sL) * sign(lidot(icross(p1,p2),closestPoint));
                        if(type == 0){//circle
                            FragColor = vec4(cL,0,sL,1);
                        }
                        if(type == 1){//horocycle
                            sL = 2*sign(sL)*mod(sL,.5);
                            FragColor = vec4(sL,0,-sL,1);
                        }
                        if(type == 2){//hypercycle
                            FragColor = vec4(2*mod(acosh(cL),.5),0,sign(sL),1);
                        }
                    } else {
                        k = (safe_sqrt(2*coshdis - 2))/(wallcornershad)/-iaunit;
                        FragColor *= vec4(k,k,k,1);
                    }
                    //...
                }

                {
                    float coshpr = cosh(pr*iaunit);
                    //...
                    float grad = lidot(endOfRay,pl);
                    if(pw == current_world && grad <= coshpr && !h2PortalBetween(pl,endOfRay,iopip)){
                        grad = -sqrt(2*grad - 2)/pr/iaunit;
                        float side = mirrorPlayer*sign(lidot(icross(pl,pl2),endOfRay));
                        if(lidot(endOfRay,pl2) < coshpr) FragColor = vec4(1,grad,side,1);
                        else FragColor = vec4(0,grad,side,1);
                        //FragColor = vec4(0,0,1,1);
                        //...
                    }
                    grad = lidot(endOfRay,playerduploc);
                    if(grad <= coshpr && duppw == current_world && h2PortalBetween(playerduploc,endOfRay,diopip)){
                        grad = -sqrt(2*grad - 2)/pr/iaunit;
                        float side = mirrorDup*sign(lidot(icross(playerduploc,playerduploc2),endOfRay));
                        if(lidot(endOfRay,playerduploc2) < coshpr) FragColor = vec4(1,grad,side,1);
                        else FragColor = vec4(0,grad,side,1);
                    }
                    //...

                }
            }
            if(iaunit == 0 ){
                //FragColor = vec4(1,0,0,1);
                vec2 closestPoint, closestClosestPoint;
                float nwdis2 = wallcornershad*wallcornershad;
                while(i < portsAndWalls[current_world + 1]){
                    int dat = int(portsAndWalls[i + 7]);
                    if((dat&4) > 0){
                        vec2 p1 = vec2(portsAndWalls[i],portsAndWalls[i+1]),
                        p2 = vec2(portsAndWalls[i+3],portsAndWalls[i+4]) - p1,
                        endOfRayCopy = vec2(endOfRay.x,endOfRay.y) - p1;
                        float rl2 = dot(p2,p2);
                        int type = dat&3;
                        if(type == 0){
                            closestPoint = endOfRayCopy*sqrt(rl2/dot(endOfRayCopy,endOfRayCopy));
                            float limit = portsAndWalls[i + 6];
                            float dis = distance(closestPoint, p2);
                            dis *= dis;
                            if(dis > limit){
                                float c = 1 - limit*0.5/rl2, s = sqrt(1 - c*c);
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
                if(distance(endOfRay,pl) < pr && current_world == pw && !e2PortalBetween(endOfRay,pl,iopip)){
                    float grad = distance(endOfRay,pl)/pr;
                    float side = mirrorPlayer*e2side(pl,pl2,endOfRay);
                    if(distance(endOfRay,pl2) < pr) FragColor=vec4(1,grad,side,1);
                    else FragColor=vec4(0,grad,side,1);
                } else if(distance(endOfRay,playerduploc) < pr && current_world == duppw && e2PortalBetween(endOfRay,playerduploc,diopip)){
                    float grad = distance(endOfRay,playerduploc)/pr;
                    float side = mirrorDup*e2side(playerduploc,playerduploc2,endOfRay);
                    if(distance(endOfRay,playerduploc2) < pr) FragColor=vec4(1,grad,side,1);
                    else FragColor=vec4(0,grad,side,1);
                } else {
                    vec2 z = vec2(endOfRay.x,endOfRay.y);
                    int N = 30, n = N;
                    while(n > 0 && dot(z,z) < 4){
                        t = z.x*z.x - z.y*z.y + endOfRay.x;
                        z.y = 2*z.x*z.y + endOfRay.y;
                        z.x = t;
                        n--;
                    }
                    float I=float(n)/float(N);
                    float test = portsAndWalls[66559];
                    //if(max(abs(endOfRay.x),abs(endOfRay.y))<tpi){
                        //FragColor = vec4((mod(endOfRay.x-tpi,2*tpi)-tpi)/6, (.5+test/2)*float((int(floor(endOfRay.x))%2+int(floor(endOfRay.y))%2+6)%2)/2+I/2+0.01,(mod(endOfRay.y-tpi,2*tpi)-tpi)/6,1);
                        FragColor = vec4((mod(endOfRay.x-tpi,2*tpi)-tpi)/6, (.5+test/2)*Mod(floor(endOfRay.x) + floor(endOfRay.y),2)/2+I/2+0.01,(mod(endOfRay.y-tpi,2*tpi)-tpi)/6,1);
                    //}else FragColor=vec4(0.5,0.5,0.5,1);
                    /*{
                        vec2 P = vec2(endOfRay.x,endOfRay.y);
                        P = mod(P + 6, 12) - 6;
                        FragColor = vec4((P.x+current_world)/6, (.5+test/2)*float((int(floor(P.x))%2+int(floor(P.y+current_world))%2+6)%2)/2+I/2+0.01,(P.y)/6,1);
                    }//*/
                    if(nearstWallIndex > 0){
                        i = nearstWallIndex;
                        vec2 p1 = vec2(portsAndWalls[i],portsAndWalls[i+1]),
                        p2 = vec2(portsAndWalls[i+3],portsAndWalls[i+4]) - p1,
                        endOfRayCopy = vec2(endOfRay.x,endOfRay.y) - p1;
                        if(nwdis2 < wallthickness*wallthickness){
                            float relLoc;
                            int type = int(portsAndWalls[i + 7])&3;
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
                float cospr = cos(pr*iaunit);
                float k=sign(endOfRay.x*endOfRay.y*endOfRay.z);
                vec3 color = (sign(endOfRay) + 1 - endOfRay)/2;
                FragColor = vec4(color.x,color.y,color.z,1);
                //
                vec3 closestPoint, closestClosestPoint;
                float cosWallCornerShad = cos(wallcornershad*iaunit);
                float codtccp = cosWallCornerShad;//cos of distance to closest closest point, I am not typing that all out outside of this comment
                while(i < portsAndWalls[current_world + 1]){
                    if( (int(portsAndWalls[i + 7])&4) > 0){
                        vec3 p1 = vec3(portsAndWalls[i],portsAndWalls[i + 1],portsAndWalls[i + 2]),
                        p2 = vec3(portsAndWalls[i + 3],portsAndWalls[i + 4],portsAndWalls[i + 5]);
                        float rc = dot(p1,p2), rss = 1 - rc*rc;//cos of r and sin(r)^2
                        float Retro = dot(endOfRay,p1),
                        A = safe_sqrt( rss/(1 - Retro*Retro ) ),
                        B = rc - A * Retro;
                        closestPoint = A*endOfRay + B*p1;
                        float limit = portsAndWalls[i + 6];
                        if(dot(closestPoint,p2) < limit){
                            vec3 ep1 = vec3(portsAndWalls[i + 8],portsAndWalls[i + 9],portsAndWalls[i + 10]);
                            vec3 ep2 = vec3(portsAndWalls[i + 11],portsAndWalls[i + 12],portsAndWalls[i + 13]);
                            if(dot(ep1,endOfRay) > dot(ep2,endOfRay)) closestPoint = ep1;
                            else closestPoint = ep2;
                        }
                        float cosDisToCP = dot(endOfRay, closestPoint);
                        if(cosDisToCP > codtccp*1.000001){//no "z" fighting and bounds checking
                            nearstWallIndex = i;
                            closestClosestPoint = closestPoint;
                            codtccp = cosDisToCP;
                        }
                    }
                    i += sizeOfpow;
                }
                float cosWallThickness = cos(wallthickness*iaunit);
                if(nearstWallIndex > 0) {
                    if(codtccp > cosWallThickness){
                        vec3 p1 = vec3(portsAndWalls[nearstWallIndex],portsAndWalls[nearstWallIndex + 1],portsAndWalls[nearstWallIndex + 2]),
                        p2 = vec3(portsAndWalls[nearstWallIndex + 3],portsAndWalls[nearstWallIndex + 4],portsAndWalls[nearstWallIndex + 5]);
                        float relLoc = dot(endOfRay,p2) * sign(dot(cross(p1,p2),endOfRay));
                        FragColor = vec4(relLoc,1,-relLoc,1);
                    } else {
                            k = (safe_sqrt( 2 - 2*codtccp))/(wallcornershad)/iaunit;
                            FragColor *= vec4(k,k,k,1);
                    }
                }
                //...
                if(dot(endOfRay,pl) > cospr && current_world==pw && !s2PortalBetween(endOfRay,pl,iopip)){
                    float grad = smallArccos(dot(endOfRay,pl))/iaunit/pr;
                    float side = -mirrorPlayer*sign(dot(cross(pl,pl2),endOfRay));
                    if(dot(endOfRay,pl2) > cospr) FragColor = vec4(1,grad,side,1);
                    else FragColor=vec4(0,grad,side,1);
                } else if(dot(endOfRay,playerduploc) > cospr && current_world==duppw && s2PortalBetween(endOfRay,playerduploc,diopip)){
                    float grad = smallArccos(dot(endOfRay,playerduploc))/pr/iaunit;
                    float side = -mirrorDup*sign(dot(cross(playerduploc,playerduploc2),endOfRay));
                    if(dot(endOfRay,playerduploc2) > cospr) FragColor = vec4(1,grad,side,1);
                    else FragColor = vec4(0,grad,side,1);
                }

            }
            //...
        }
    }
    if(debug) FragColor = dcolor;
    //...
}
