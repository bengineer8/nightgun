//to compile for windows on Linux: x86_64-w64-mingw32-g++ -O3 -w code.cpp -o nightgun.exe -static -L. -l glfw3dll
//to compile on Windows: g++ -O3 -w code.cpp -o nightgun.exe -static -L. -l glfw3dll
//^to get g++ working on Windows, follow this but only inclusively to "Check your MinGW installation": https://code.visualstudio.com/docs/cpp/config-mingw
//to run, open the folder this file is in, alt+d, type "cmd", hit enter, copy the compile command, and hit enter

//for linux: g++ -O3 -w code.cpp -o nightgun -lglfw
#include <iostream>
#include "glad.h"
#include "glad.c"
#include "glfw3.h"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

//#include <windows.h>

//al = r^2 (2-2cos(a)) for E2 and S2 and cosh(r)^2 - sinh(r)^2 cos(a) for H2

#define positionsSize 1048576

//SDL_Joystick* gGameController = NULL;

struct entity{
    double pos1[3];
    double pos2[3];
    int props1;
    int props2;
};


/*std::vector<float> etfl(entity E){
    std::vector<float> F={E.pos1[0],E.pos1[1],E.pos1[2],E.pos2[0],E.pos2[1],E.pos2[2], E.props1, E.props2};
    return F;
}*/

std::vector<std::vector<std::vector<float>>> pawbuffer;
//std::vector<std::vector<std::vector<double>>> entitybuffer;//things that move around are stored with more precision
//std::vector<std::vector<std::vector<double>>> dopgangbuf;//doppleganger buffer, used for render and collision checking copies, not used yet
std::vector<float> worldCurvatures;
float shader_data[positionsSize];
float worldCurvs[16];//to do: make this all done in 1 batch later

const float pi=3.14159265359, sr2=sqrt(2), isr2=sr2/2,ps=0.032,pr=1.0/9.0, ped=1.0/1000, inf = 1/0.0;//player speed. player radius, portal ejection distance
const double playerWallCollisionr = 1.0/18 + pr;
const double cosplayerWallCollisionr = cos(playerWallCollisionr);
const double coshplayerWallCollisionr = cosh(playerWallCollisionr);
const float prpwt=1.0/18+pr, prppwtc=cos(prpwt), prppwts=sin(prpwt);
const double prc=cos(pr), prs=sin(pr), pedc=cos(ped), peds=sin(ped);
const double prch = cosh(pr), prsh = sinh(pr);
const int sizeOfPow = 14;
const char *vertexShaderSource = "#version 430 core\n"
    "layout (location = 0) in vec3 aPos;\n"
    "void main()\n"
    "{\n"
    "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
    "}\0";
float plo[3];//player orientation
float pl[3],camRef[3],plp[1],duppl[3],paw[512];//player location, ..., player properties(mirror view flag (flip y)), duplicate player location, used for display purposes ,portals and walls
float pl2[3],duppl2[3];
float mirrorDup = 1;
float fortyfivedegrot[][3]={{isr2,isr2,0},{-isr2,isr2,0},{0,0,1}};
int pw, duppw, iopip = -1, diopip = -1;//index of portal intersecting player

void printM(double m[][3]){
    for(int n=0;n<3;n++){
        printf("%.9lf\t%.9lf\t%.9lf\n",m[0][n],m[1][n],m[2][n]);
    }
    printf("\n");
}

void printP(double p[3]){
    printf("%.9lf, %.9lf, %.9lf\n",p[0],p[1],p[2]);
}

void printP(float p[3]){
    printf("%.9g, %.9g, %.9g\n",p[0],p[1],p[2]);
}

float mkdest(int t, int w, int i, int s, int m) {return (float)(t|((w|((i|((s|m<<1)<<10))<<9))<<3));}//type, destination (world, index, side, mirror)

double safe_sqrt(double x){
    if(x<0) return 0;
    return sqrt(x);
}

char sign(float f){
    if(f<0) return -1;
    return 1;
}

//comment this out on mac
/*double abs(double x){
    if(x<0) return -x;
    return x;
}//*/

char sign(double f){
    if(f<0) return -1;
    return 1;
}

double arctan(double s, double c){
    double t=-atan(c/s);
    if(s<0) c=1.5;
    else c=0.5;
    return pi*c+t;
}


double disSquaredXY(double a[], double b[]){
    double c[2] = {a[0] - b[0], a[1] - b[1]};
    return c[0]*c[0] + c[1]*c[1];
}

void vec3(double v[], double x, double y, double z){
    v[0]=x;v[1]=y;v[2]=z;
}

void vec3(float v[], double x, double y, double z){
    v[0]=x;v[1]=y;v[2]=z;
}

void copypt(double p1[3], double p2[3]){
    p2[0]=p1[0];p2[1]=p1[1];p2[2]=p1[2];
}

void copypt(float p1[3], double p2[3]){
    p2[0]=p1[0];p2[1]=p1[1];p2[2]=p1[2];
}

void copypt(double p1[3], float p2[3]){
    p2[0]=p1[0];p2[1]=p1[1];p2[2]=p1[2];
}

void copymat(double a[3][3], double b[3][3]){
    copypt(a[0],b[0]);
    copypt(a[1],b[1]);
    copypt(a[2],b[2]);
}

void mat3(double m[3][3], double a, double b, double c, double d, double e, double f, double g, double h, double i){//same format as in GLSL for ease of porting
    m[0][0]=a;m[0][1]=b;m[0][2]=c;
    m[1][0]=d;m[1][1]=e;m[1][2]=f;
    m[2][0]=g;m[2][1]=h;m[2][2]=i;
}

void matxpt(double m[][3], double p[3]){
    double A=m[0][0]*p[0]+m[1][0]*p[1]+m[2][0]*p[2],
    B=m[0][1]*p[0]+m[1][1]*p[1]+m[2][1]*p[2];
    p[2]=m[0][2]*p[0]+m[1][2]*p[1]+m[2][2]*p[2];
    p[0]=A;
    p[1]=B;
}

void matxpt(double m[][3], float p[3]){
    double A=m[0][0]*p[0]+m[1][0]*p[1]+m[2][0]*p[2],
    B=m[0][1]*p[0]+m[1][1]*p[1]+m[2][1]*p[2];
    p[2]=m[0][2]*p[0]+m[1][2]*p[1]+m[2][2]*p[2];
    p[0]=A;
    p[1]=B;
}

void matxpt(double m[][3], double p[3], double p2[3]){
    p2[0]=m[0][0]*p[0]+m[1][0]*p[1]+m[2][0]*p[2];
    p2[1]=m[0][1]*p[0]+m[1][1]*p[1]+m[2][1]*p[2];
    p2[2]=m[0][2]*p[0]+m[1][2]*p[1]+m[2][2]*p[2];
}

void matxpt(double m[][3], float p[3], float p2[3]){
    p2[0]=m[0][0]*p[0]+m[1][0]*p[1]+m[2][0]*p[2];
    p2[1]=m[0][1]*p[0]+m[1][1]*p[1]+m[2][1]*p[2];
    p2[2]=m[0][2]*p[0]+m[1][2]*p[1]+m[2][2]*p[2];
}

void transpose(double a[][3], double b[][3]){//inflicts the transpose of a on to b
    b[0][0]=a[0][0];b[1][1]=a[1][1];b[2][2]=a[2][2];
    b[0][1]=a[1][0];b[0][2]=a[2][0];b[1][2]=a[2][1];
    b[1][0]=a[0][1];b[2][0]=a[0][2];b[2][1]=a[1][2];
}//this is also how you invert rotation+translation in s2

void rotate(double& x, double& y, double c, double s){
    double t = x*c - y*s;
    y = x*s + y*c;
    x=t;
}

void rotate(float& x, float& y, double c, double s){
    double t = x*c - y*s;
    y = x*s + y*c;
    x=t;
}

void rotXY(double m[][3], double c, double s){//to do, swap order of inputs
    rotate(m[0][0],m[0][1],c,s);
    rotate(m[1][0],m[1][1],c,s);
    rotate(m[2][0],m[2][1],c,s);
}

void rotXZ(double m[][3], double c, double s){
    rotate(m[0][0],m[0][2],c,s);
    rotate(m[1][0],m[1][2],c,s);
    rotate(m[2][0],m[2][2],c,s);
}

void s2matto(double p[3], double m[][3]){//turns m into a transformation matrix to translate p to (0,0,1) on a sphere
    double r=sqrt(p[0]*p[0]+p[1]*p[1]), c=p[0]/r, s=p[1]/r;
    if(r==0) {c=1;s=0;}
    mat3(m, p[2]*c,-s,p[0], p[2]*s,c,p[1], -r,0,p[2]);
}

void s2matto(float p[3], double m[][3]){//turns m into a transformation matrix to translate p to (0,0,1) on a sphere
    double r=sqrt(p[0]*p[0]+p[1]*p[1]), c=p[0]/r, s=p[1]/r;
    if(r==0) {c=1;s=0;}
    mat3(m, p[2]*c,-s,p[0], p[2]*s,c,p[1], -r,0,p[2]);
}

void s2matto(int i, double m[][3]){
    double x = paw[i],
    y = paw[i + 1],
    z = paw[i + 2],
    ns = paw[i + 8],
    c = paw[i + 9],
    nr = paw[i + 10],
    c2 = paw[i + 11],
    s2 = paw[i + 12];
    mat3(m, z*c,ns,x, -z*ns,c,y, nr,0,z);
    rotXY(m,c2,-s2);
}

void s2dadtopt(double d, double c, double s, double p[3]){//distance and direction to 3D pt
    p[0]=sin(d);
    p[1]=p[0]*s;
    p[0]*=c;
    p[2]=cos(d);
}

void s2dadtopt(double d, double c, double s, float p[3]){//distance and direction to 3D pt
    p[0]=sin(d);
    p[1]=p[0]*s;
    p[0]*=c;
    p[2]=cos(d);
}

double s2disrank(double c,double s){
    c+=1;
    if(s<0) c=-1;
    return 2-c;
}

void normFix2D(float v[]){//makes a "unit" vector more of a "unit" vector
    float r=sqrt(v[0]*v[0]+v[1]*v[1]);
    v[0]/=r;v[1]/=r;
    if(v[2]>1) v[2]=sr2;
}

void backOnSphere(float p[]){
    float r=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[0]/=r;p[1]/=r;p[2]/=r;
}

void backOnSphere(double p[]){
    float r=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[0]/=r;p[1]/=r;p[2]/=r;
}

double dot(double p1[3], double p2[3]){
    return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
}

double dot(float p1[3], double p2[3]){
    return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
}

void cross(double p1[3], double p2[3], double p3[3]){
    p3[0] = p1[1]*p2[2] - p1[2]*p2[1];
    p3[1] = p1[2]*p2[0] - p1[0]*p2[2];
    p3[2] = p1[0]*p2[1] - p1[1]*p2[0];
}

double s2side(double p1[3], double p2[3], double p3[3]){
    double C[3];
    cross(p1,p2,C);
    if(dot(p3,C) < 0) return -1;
    else return 1;
}

double lidot(double a[], double b[]){
    return -a[0]*b[0] - a[1]*b[1] + a[2]*b[2];
}

void h2invert(double a[][3], double b[][3]){//makes b the inverse of a if a is an isometry of H2    ...that sounds way too technical for what it does
    b[0][0]=a[0][0];b[1][1]=a[1][1];b[2][2]=a[2][2];
    b[0][1]=a[1][0];b[0][2]=-a[2][0];b[1][2]=-a[2][1];
    b[1][0]=a[0][1];b[2][0]=-a[0][2];b[2][1]=-a[1][2];
}//this is also how you invert rotation+translation in h2

void h2matto(double p[], double m[][3]){
    double r=sqrt(p[0]*p[0]+p[1]*p[1]), c=p[0]/r, s=p[1]/r;
    if(r==0){
        c=1;
        s=0;
    }
    mat3(m, p[2]*c,-s,-p[0], p[2]*s,c,-p[1], -r,0,p[2]);
}

void h2matto(float p[], double m[][3]){
    double r=sqrt(p[0]*p[0]+p[1]*p[1]), c=p[0]/r, s=p[1]/r;
    if(r==0){
        c=1;
        s=0;
    }
    mat3(m, p[2]*c,-s,-p[0], p[2]*s,c,-p[1], -r,0,p[2]);
}

void lorenz(double& x, double& y, double c, double s){
    double t = x*c + y*s;
    y = x*s + y*c;
    x=t;
}

void lorenz(float& x, float& y, double c, double s){
    double t = x*c + y*s;
    y = x*s + y*c;
    x=t;
}

void lorenzXZ(double m[][3], double c, double s){
    lorenz(m[0][0],m[0][2],c,s);
    lorenz(m[1][0],m[1][2],c,s);
    lorenz(m[2][0],m[2][2],c,s);
}

void lorenzYZ(double m[][3], double c, double s){
    lorenz(m[0][1],m[0][2],c,s);
    lorenz(m[1][1],m[1][2],c,s);
    lorenz(m[2][1],m[2][2],c,s);
}

void idealrot(double m[3][3], double d){
    double d1 = d*d/2, d2 = 1 - d1, d3 = 1 + d1, a, b, c;
    for(int n = 0; n < 3; n++){
        a = m[n][0]*d2 - m[n][1]*d + m[n][2]*d1;
        b = (m[n][0] - m[n][2])*d + m[n][1];
        c = -m[n][0]*d1 - m[n][1]*d + m[n][2]*d3;
        m[n][0] = a; m[n][1] = b; m[n][2] = c;
    }
}

void h2dadtopt(double d, double c, double s, double p[3]){//distance and direction to 3D pt
    p[0]=sinh(d);
    p[1]=p[0]*s;
    p[0]*=c;
    p[2]=cosh(d);
}

void h2dadtopt(double d, double c, double s, float p[3]){//distance and direction to 3D pt
    p[0]=sinh(d);
    p[1]=p[0]*s;
    p[0]*=c;
    p[2]=cosh(d);
}

void backOnHyperboloid(float p[]){
    p[2] = sqrt(p[0]*p[0] + p[1]*p[1] + 1);
}

double h2side(double p1[3], double p2[3], double p3[3]){
    double icross[3] = {p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2], p1[1]*p2[0] - p1[0]*p2[1]};
    if(lidot(icross,p3) < 0) return -1;
    return 1;
}

void licross(double p1[3], double p2[3], double p3[3]){
    vec3(p3, p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2], p1[1]*p2[0] - p1[0]*p2[1]);
}

bool h2PortalBetween(double A[3], double B[3], int i){
    if(i > 0){
        double p1[3] = {paw[i], paw[i + 1], paw[i + 2]}, p2[3] = {paw[i + 3], paw[i + 4], paw[i + 5]};
        double k = lidot(p1,p2), iA = lidot(A,p1), iB = lidot(B,p1), s = k - iB;
        k *= s;
        if(iA*s > k && k > s*iB) return true;
    }
    return false;
}

void addExtraToPOW2(std::vector<float>& G, float iaunit){
    double p1[3] = {G[0],G[1],G[2]}, p2[3] = {G[3],G[4],G[5]};
    double L = G[6];
    int dat = (int)G[7];
    char type = dat&7;
    bool iswall = (type > 3);
    if(iaunit > 0){
        if(iswall){
            double rot[3][3], roti[3][3];
            double endpt1[3], endpt2[3];
            s2matto(p1,rot);
            matxpt(rot,p2);
            double cos = G[6] + 1, sin = sqrt(1 - cos*cos);
            copypt(p2,endpt1);
            copypt(p2,endpt2);
            rotate(endpt1[0],endpt1[1],cos,sin);
            rotate(endpt2[0],endpt2[1],cos,-sin);
            transpose(rot,roti);
            matxpt(roti,endpt1);
            matxpt(roti,endpt2);
            G.push_back(endpt1[0]);
            G.push_back(endpt1[1]);
            G.push_back(endpt1[2]);
            G.push_back(endpt2[0]);
            G.push_back(endpt2[1]);
            G.push_back(endpt2[2]);
        } else {
            double rc = dot(p1,p2),
            rs = sqrt(1 - rc*rc);
            double rot[3][3], roti[3][3];
            s2matto(p1,rot);
            matxpt(rot,p2);
            G.push_back(rot[0][1]);//-s1
            G.push_back(rot[1][1]);//c1
            G.push_back(rot[2][0]);//-r
            G.push_back(p2[0]/rs);//c2
            G.push_back(p2[1]/rs);//s2
            G.push_back(rs);
        }
    } else if(iaunit == 0){
            p2[0] -= p1[0];p2[1] -= p1[1];
            if(iswall){
                if(type == 5) {for(int n = 0; n < 6; n++) G.push_back(0);}
                else {
                    double cos = 1 - G[6]/2, sin = sqrt(1 - cos*cos);
                    double x = p2[0], y = p2[1];
                    rotate(x,y,cos,sin);
                    G.push_back(x);G.push_back(y);G.push_back(0);
                    x = p2[0]; y = p2[1];
                    rotate(x,y,cos,-sin);
                    G.push_back(x);G.push_back(y);G.push_back(0);
                }
        } else {
            //Note: unused ATM
            p2[0] -= p1[0];p2[1] -= p1[1];
            double r = sqrt(p2[0]*p2[0] + p2[1]*p2[1]);
            G.push_back(p2[0]/r);
            G.push_back(-p2[1]/r);
            G.push_back(r);
            double ox = -p1[0], oy = -p1[1];
            rotate(ox,oy,G[8],G[9]);
            G.push_back(ox);
            G.push_back(oy);
            G.push_back(0);//to do: put something funny here?
        }
    } else if(iaunit < 0){
        if(iswall){
            double k = lidot(p1,p2);
            double rot[3][3], roti[3][3];
            h2matto(p2,rot);
            matxpt(rot,p1);
            double r = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
            rotXY(rot, p1[0]/r,-p1[1]/r);
            double endpt1[3], endpt2[3];
            endpt1[0] = k*(L - 1)/r;
            endpt1[1] = sqrt(L*L - 1 - endpt1[0]*endpt1[0]);
            endpt1[2] = sqrt(endpt1[0]*endpt1[0] + endpt1[1]*endpt1[1] + 1);
            copypt(endpt1,endpt2);
            endpt2[1] = -endpt2[1];
            h2invert(rot,roti);
            matxpt(roti,endpt1);
            matxpt(roti,endpt2);
            G.push_back(endpt1[0]);
            G.push_back(endpt1[1]);
            G.push_back(endpt1[2]);
            G.push_back(endpt2[0]);
            G.push_back(endpt2[1]);
            G.push_back(endpt2[2]);
            //...
        } else {
            for(int n = 8; n < 14; n++) G.push_back(0);//placeholder
            //...
        }
    }
    //...
}
/*void addExtraToPOW(std::vector<float>& G, float iaunit){//NOTE: to do: split this up by geometry instead of by type
    double p1[3] = {G[0],G[1],G[2]}, p2[3] = {G[3],G[4],G[5]};
    int dat = (int)G[7];
    char type = dat&7;
    if(type > 3){
        double endpt1[3], endpt2[3];
        if(iaunit == 0){
            if(type == 5){
                copypt(p1,endpt1);
                copypt(p2,endpt2);
            } else if(iaunit > 0){
                p2[0] -= p1[0];p2[1] -= p1[1];
                double cos = 1 - G[6]/2, sin = sqrt(1 - cos*cos);
                copypt(p2,endpt1);
                copypt(p2,endpt2);
                rotate(endpt1[0],endpt1[1],cos,sin);
                rotate(endpt2[0],endpt2[1],cos,-sin);
                endpt1[0] += p1[0];endpt1[1] += p1[1];
                endpt2[0] += p1[0];endpt2[1] += p1[1];
            }
        } else if(iaunit > 0){
            double rot[3][3], roti[3][3];
            s2matto(p1,rot);
            matxpt(rot,p2);
            double cos = G[6] + 1, sin = sqrt(1 - cos*cos);
            copypt(p2,endpt1);
            copypt(p2,endpt2);
            rotate(endpt1[0],endpt1[1],cos,sin);
            rotate(endpt2[0],endpt2[1],cos,-sin);
            transpose(rot,roti);
            matxpt(roti,endpt1);
            matxpt(roti,endpt2);
        }
        G.push_back(endpt1[0]);
        G.push_back(endpt1[1]);
        G.push_back(endpt1[2]);
        G.push_back(endpt2[0]);
        G.push_back(endpt2[1]);
        G.push_back(endpt2[2]);
    } else {
        if(iaunit == 0){
            p2[0] -= p1[0];p2[1] -= p1[1];
            double r = sqrt(p2[0]*p2[0] + p2[1]*p2[1]);
            G.push_back(p2[0]/r);
            G.push_back(-p2[1]/r);
            G.push_back(r);
            double ox = -p1[0], oy = -p1[1];
            rotate(ox,oy,G[8],G[9]);
            G.push_back(ox);
            G.push_back(oy);
            G.push_back(0);//to do: put something funny here?
        } else if(iaunit > 0){
            double rc = dot(p1,p2),
            rs = sqrt(1 - rc*rc);
            double rot[3][3], roti[3][3];
            s2matto(p1,rot);
            matxpt(rot,p2);
            G.push_back(rot[0][1]);//-s1
            G.push_back(rot[1][1]);//c1
            G.push_back(rot[2][0]);//-r
            G.push_back(p2[0]/rs);//c2
            G.push_back(p2[1]/rs);//s2
            G.push_back(rs);
        }
    }
    if(iaunit < 0){
        for(int n = 8;n < 14;n++) G.push_back(0);//placeholder
    }
}*/

void updateDuplicateRot(float pos[], float pos2[], float dpos[], float dpos2[], int w, int i){
    double p1[] = {paw[i],paw[i + 1],paw[i + 2]};
    int dat = (int)paw[i + 7];
    int type = dat&3;
    char mirror = (dat>>23)&1, side = (dat>>22)&1;
    double iaunit = worldCurvatures[w];
    double pos2r[3];
    copypt(pos2,pos2r);
    if(iaunit == 0 && type == 1){//lines in E2 are special needs
        double p0[] = {pos[0],pos[1]}, p2[] = {paw[i + 3],paw[i + 4]};
        p0[0] -= p1[0]; p0[1] -= p1[1];
        p2[0] -= p1[0]; p2[1] -= p1[1];
        double C = (p0[0]*p2[0] + p0[1]*p2[1])/(p2[0]*p2[0] + p2[1]*p2[1]);
        p2[0] *= C; p2[1] *= C;
        p1[0] += p2[0]; p1[1] += p2[1];
    }
    if(iaunit == 0){
        p1[0] -= pos[0];p1[1] -= pos[1];
        pos2r[0] -= pos[0];pos2r[1] -= pos[1];
    } else if(iaunit > 0){
        double rot[3][3];
        s2matto(pos,rot);
        matxpt(rot,p1);
        matxpt(rot,pos2r);
    } else if(iaunit < 0){
        double rot[3][3];
        h2matto(pos,rot);
        matxpt(rot,p1);
        matxpt(rot,pos2r);
    }
    p1[2] = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
    p1[0] /= p1[2]; p1[1] /= p1[2];
    rotate(pos2r[0],pos2r[1],p1[0],-p1[1]);
    if(mirror == 1) pos2r[1] = -pos2r[1];
    if(side == 1){
        pos2r[1] = -pos2r[1];
        pos2r[0] = -pos2r[0];
    }
    if(iaunit == 0){
        pos2r[0] /= pr; pos2r[1] /= pr;
    } else if(iaunit > 0){
        pos2r[0] /= prs; pos2r[1] /= prs;
    } else if(iaunit < 0){
        pos2r[0] /= prsh; pos2r[1] /= prsh;
    }
    w = (dat>>3)&511;
    i = int(paw[w]) + sizeOfPow*(dat>>12)&1023;
    iaunit = worldCurvatures[w];
    vec3(p1,paw[i],paw[i + 1],paw[i + 2]);
    if(iaunit == 0 && type == 1){
        double p0[] = {dpos[0],dpos[1]}, p2[] = {paw[i + 3],paw[i + 4]};
        p0[0] -= p1[0]; p0[1] -= p1[1];
        p2[0] -= p1[0]; p2[1] -= p1[1];
        double C = (p0[0]*p2[0] + p0[1]*p2[1])/(p2[0]*p2[0] + p2[1]*p2[1]);
        p2[0] *= C; p2[1] *= C;
        p1[0] += p2[0]; p1[1] += p2[1];
    }
    if(iaunit == 0){
        pos2r[0] *= pr; pos2r[1] *= pr;
        p1[0] -= dpos[0];p1[1] -= dpos[1];
        p1[2] = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
        p1[0] /= p1[2]; p1[1] /= p1[2];
        rotate(pos2r[0],pos2r[1],p1[0],p1[1]);
        pos2r[0] += dpos[0]; pos2r[1] += dpos[1];
        dpos2[0] = pos2r[0]; dpos2[1] = pos2r[1];
    } else if(iaunit > 0){
        double rot[3][3], roti[3][3];
        s2matto(dpos,rot);
        matxpt(rot,p1);
        p1[2] = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
        p1[0] /= p1[2]; p1[1] /= p1[2];
        s2dadtopt(pr,pos2r[0],pos2r[1],dpos2);
        rotate(dpos2[0],dpos2[1],p1[0],p1[1]);
        transpose(rot,roti);
        matxpt(roti,dpos2);
    } else if(iaunit < 0){
        double rot[3][3], roti[3][3];
        h2matto(dpos,rot);
        matxpt(rot,p1);
        p1[2] = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
        p1[0] /= p1[2]; p1[1] /= p1[2];
        h2dadtopt(pr,pos2r[0],pos2r[1],dpos2);
        rotate(dpos2[0],dpos2[1],p1[0],p1[1]);
        h2invert(rot,roti);
        matxpt(roti,dpos2);
        //...
    }
    //
}


void createDuplicate(double RL[], double RD, double p[], int ipi, int* w, double pos1[], double pos2[]){
    //RL[2] = RD;
    //printP(RL);
    //
    mirrorDup = plp[0];//TODO make this entity specific
    int dat = (int)paw[ipi + 7];
    int type = dat & 3;
    char mirror = (dat>>23)&1, side=(dat>>22)&1;
    if(mirror==1) {//mirror
        RL[1]=-RL[1];
        p[1]=-p[1];
        mirrorDup = -mirrorDup;//TODO make this entity specific
    }
    if(side==1) {//side swap
        RL[1]=-RL[1];
        RD=-RD;
        p[0]=-p[0];
        p[1]=-p[1];
    }
    *w=((dat)>>3)&511;
    double iaunit=worldCurvatures[*w];
    int i=int(paw[*w])+sizeOfPow*(dat>>12)&1023;
    diopip = i;//TODO make this entity specific
    double p1[]={paw[i],paw[i+1],paw[i+2]},p2[]={paw[i+3],paw[i+4],paw[i+5]};
    if(iaunit == 0){
        pos1[2]=0;
        pos2[2]=0;
        p2[0]-=p1[0];p2[1]-=p1[1];
        double r=sqrt(p2[0]*p2[0]+p2[1]*p2[1]);
        p2[0]/=r;p2[1]/=r;
        switch(type){
            case 0:
                pos1[0]=r+RD;
                pos1[1]=pos1[0]*RL[1];
                pos1[0]=pos1[0]*RL[0];
                rotate(p[0],p[1],RL[0],RL[1]);
                break;
            case 1:
                {
                    double t=p[0];
                    p[0]=-p[1];
                    p[1]=t;
                }
                pos1[0]=r/2 - RL[1];
                pos1[1]=RD;
                break;
            default:
                break;
        }
        pos2[0]=p[0]*p[2]+pos1[0];pos2[1]=p[1]*p[2]+pos1[1];
        rotate(pos1[0],pos1[1],p2[0],p2[1]);
        rotate(pos2[0],pos2[1],p2[0],p2[1]);
        pos1[0]+=p1[0];pos1[1]+=p1[1];
        pos2[0]+=p1[0];pos2[1]+=p1[1];
    }
    if(iaunit > 0){
        double rc=dot(p1,p2), rs=safe_sqrt(1-rc*rc),c=cos(RD),s=sin(RD);
        double rot[3][3], roti[3][3];
        s2matto(i,rot);
        vec3(pos1, rs,0,rc);
        vec3(pos2, prs*p[0], prs*p[1], prc);
        rotate(pos2[2],pos2[0],rc,rs);
        rotate(pos2[2],pos2[0],c,s);
        rotate(pos1[2],pos1[0],c,s);
        rotate(pos1[0],pos1[1],RL[0],RL[1]);
        rotate(pos2[0],pos2[1],RL[0],RL[1]);
        transpose(rot,roti);
        matxpt(roti,pos1);
        matxpt(roti,pos2);
    }
    if(iaunit < 0){
        if(type > 0) {//bandaid TODO actually fix this
            RD = -RD;
            p[0] = -p[0];
            p[1] = -p[1];
        }
        vec3(pos1, sinh(RD),0,cosh(RD));
        vec3(pos2, prsh*p[0],prsh*p[1],prch);
        lorenz(pos2[0],pos2[2], pos1[2],pos1[0]);
        //...
        double k = lidot(p1,p2);
        double rot[3][3];
        if(type > 0){
            double tv[3];
            copypt(p1,tv);
            copypt(p2,p1);
            copypt(tv,p2);
        }
        h2matto(p1,rot);
        matxpt(rot,p2);
        double r = sqrt(p2[0]*p2[0] + p2[1]*p2[1]);
        p2[0] /= r; p2[1] /= r;
        rotXY(rot, p2[0],-p2[1]);
        switch(type){
            case 0:
                rotXY(rot, RL[0],-RL[1]);
                lorenzXZ(rot,k,-r);
                break;
            case 1:
                idealrot(rot,-RL[1]);
                break;
            case 2:
                lorenzXZ(rot,r,-k);
                lorenzYZ(rot,RL[0],RL[1]);
                lorenzXZ(rot,r,k);
                break;
            default:
                printf("invalid type\n");
        }
        double roti[3][3];
        h2invert(rot,roti);
        matxpt(roti,pos1);
        matxpt(roti,pos2);
    }
}

void moveEntity(double dx, double dy, float pos[], float pos2[], float ref[], int* cw, float props[]){
    int maxits = 19;
    bool bounce = 0;//temp name, disables wall colliders and portal dupe creation because they use the same code. Makes you bounce off walls because that is what they do as a failsafe for if you go too fast to get past their colliders.
    //printf("dx:%.19lf \t dy:%.19lf \t pos:%.19f, %.19f, %.19f \t cam:%.19f, %.19f, %.19f \t cw:%i\n", dx, dy, pos[0], pos[1], pos[2], ref[0], ref[1], ref[2], *cw);
    //WIP
    //TO DO: Make it so all internal calculations are doubles
    double ogpos[3], ogpos2[3], ogref[3];
    int ogcw = *cw;
    copypt(pos,ogpos);
    copypt(pos2,ogpos2);
    copypt(ref,ogref);
    duppw = -1;//temporary hard coding, make it only do this for the entity specific duplicat later
    iopip = -1;
    diopip = -1;
    char reflect=0;//temporary
    char type = 0,mirror=1,side=0;
    int si=-1;
    if(props[0]==-1) {dy=-dy;}//correct the confusion of the player
    double d=sqrt(dx*dx+dy*dy);
    dx/=d;dy/=d;
    double sp1[]={pos[0],pos[1],pos[2]}, sp2[]={ref[0],ref[1],ref[2]};
    double EL[2]={dx,-dy}, EA[2];//Entrance/exit location and angle
    double pos2r[3];
    int iterations=0;
    double rotmstart[3][3];
    bool debug = false;
    while(d>0 && iterations < maxits){
        float iaunit = worldCurvatures[*cw];
        d+=ped;
        iterations++;
        //if(iterations>5) abort();
        if(iaunit==0){//e2
            double off[2], rot[2];
            double r,t;
            if(type < 4){
                off[0] = -sp1[0];off[1] = -sp1[1];
                sp2[0]-=sp1[0];sp2[1]-=sp1[1];
                r=sqrt(sp2[0]*sp2[0]+sp2[1]*sp2[1]);
                rot[0] = sp2[0]/r; rot[1] = -sp2[1]/r;//rotation vector to bring sp2 to (1,0)
                if(type==0) rotate(rot[0],rot[1],EL[0],EL[1]);
                rotate(off[0],off[1],rot[0],rot[1]);//apply rotations to what will be the offset
                if(si > 0){
                    if(type==0) off[0]-=r;
                    if(type==1) {//E2 lines are special needs
                        off[0]-=EL[1]+r/2;
                        t=EA[0];
                        EA[0]=-EA[1];
                        EA[1]=t;
                    }
                    rotate(off[0],off[1],EA[0],EA[1]);
                    rotate(rot[0],rot[1],EA[0],EA[1]);
                } else {
                    pos2r[0]=pos2[0];pos2r[1]=pos2[1];
                    rotate(pos2r[0],pos2r[1],rot[0],rot[1]);
                    pos2r[0]+=off[0];pos2r[1]+=off[1];
                    pos2r[2]=sqrt(pos2r[0]*pos2r[0]+pos2r[1]*pos2r[1]);
                    pos2r[0]=pos2r[0]/pos2r[2];pos2r[1]=pos2r[1]/pos2r[2];
                    pos2r[2]=pr;
                }
            } else {
                rot[0] = rotmstart[0][0]; rot[1] = rotmstart[0][1];
                off[0] = rotmstart[1][0]; off[1] = rotmstart[1][1];
            }
            int i=(int)paw[*cw],cii=-1;
            double cip1[]={0,0},cip2[]={0,0};
            double cid=d;
            //portal collision code below to do: add walls to this
            while(i<paw[*cw+1]){
                double p1[]={paw[i],paw[i+1]},p2[]={paw[i+3],paw[i+4]};
                r=(p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]);//actually r^2
                double al=paw[i+6];
                rotate(p1[0],p1[1],rot[0],rot[1]);
                p1[0]+=off[0];p1[1]+=off[1];
                rotate(p2[0],p2[1],rot[0],rot[1]);
                p2[0]+=off[0];p2[1]+=off[1];
                //if(i==10) printf("%lf, %lf, %lf, %lf\n", p1[0], p1[1], p2[0], p2[1]);
                type = ((char)paw[i+7])&3;
                if(type == 0){
                    if(p1[1]*p1[1]<=r&&(al==4*r||p2[1]*p2[1]<=al*r)){
                        double o=sqrt(r-p1[1]*p1[1]), xt=p1[0]-o, xt1=abs(xt);
                        if(((al==4) || ((xt-p2[0])*(xt-p2[0])+p2[1]*p2[1]<=al*r)) && xt>0 && xt<cid && i!=si ){
                            cip1[0]=p1[0];cip1[1]=p1[1];
                            cip2[0]=p2[0];cip2[1]=p2[1];
                            cii=i;cid=xt;
                        }
                        xt=p1[0]+o;
                        //printf("%.99lf\t%lf\n",xt,xt1);
                        if(xt>0){
                            if(((al==4) || ((xt-p2[0])*(xt-p2[0])+p2[1]*p2[1]<=al*r)) && xt<cid && (i!=si || xt>xt1) ){//problem place, think it is fixed now.
                                //std::cout<<"\t"<<((xt-p1[0])*(xt-p1[0])+p1[1]*p1[1]+ped-r)<<"\n";
                                cip1[0]=p1[0];cip1[1]=p1[1];
                                cip2[0]=p2[0];cip2[1]=p2[1];
                                cii=i;cid=xt;
                            }
                        }
                    }
                } else {
                    if((p1[1]*p2[1]<0)&&((p1[0]>=0)||(p2[0]>=0))&&(i!=si)){//crosses x axis and one of the 2 points is in front and not the source portal
                        double xt=(p1[0]-p2[0])*p2[1]/(p2[1]-p1[1])+p2[0];
                        if(xt>0 && xt<cid){
                            cip1[0]=p1[0];cip1[1]=p1[1];
                            cip2[0]=p2[0];cip2[1]=p2[1];
                            cii=i;cid=xt;
                        }
                    }
                }
                i+=sizeOfPow;
            }
            //if(cii==si) std::cout<<cid<<"\n";
            if(cii>-1){
                int dat = (int)paw[cii + 7];
                type = dat&7;
                if(type > 3){
                    si = cii;
                    d -= cid;
                    off[0] -= cid;
                    cip1[0] -= cid;
                    if(type == 5){
                        t = cip1[0];
                        cip1[0] = -cip1[1];
                        cip1[1] = t;
                    }
                    double r2 = cip1[0]*cip1[0] + cip1[1]*cip1[1];
                    double c = 2*cip1[1]*cip1[1]/r2 - 1, s = -2*cip1[0]*cip1[1]/r2;
                    rotate(rot[0],rot[1],c,-s);
                    rotate(off[0],off[1],c,-s);
                    rotate(pos2r[0],pos2r[1],c,-s);
                    rotmstart[0][0] = rot[0]; rotmstart[0][1] = rot[1];
                    rotmstart[1][0] = off[0]; rotmstart[1][1] = off[1];
                } else {
                    //handle portal collisions
                    if(((dat>>22)&1)==1) side=-1;
                    else side=1;
                    if(((dat>>23)&1)==1) {mirror=-1;props[0]=-props[0];}//might have issues
                    else mirror=1;
                    cip1[0]-=cid;
                    r=sqrt(cip1[0]*cip1[0]+cip1[1]*cip1[1]);
                    t=cip1[0]/r*side;
                    if(type==1) t*=mirror;
                    EA[0]=-t;
                    t=cip1[1]/r*side;
                    if(type==0) t*=mirror;
                    EA[1]=-t;
                    if(reflect){
                        EA[0]=-EA[0];
                        EA[1]=-EA[1];
                    }
                    if(type==0){
                        cip2[0]-=cid;
                        cip2[0]=(cip2[0]-cip1[0]);
                        cip2[1]=-(cip2[1]-cip1[1]);
                        cip1[0]*=-1;cip1[1]*=-1;
                        rotate(cip1[0],cip1[1],cip2[0],cip2[1]);
                        r=sqrt(cip1[0]*cip1[0]+cip1[1]*cip1[1]);
                        cip1[0]/=r;
                        cip1[1]/=r*mirror*side;
                        EL[0]=cip1[0];EL[1]=-cip1[1];
                    }
                    if(type==1){
                        cip1[0]-=cip2[0];cip1[1]-=cip2[1];
                        EL[1]=(r-sqrt(cip1[0]*cip1[0]+cip1[1]*cip1[1])/2)*mirror*side;
                        t=-EA[0];//this is for when it will need to link to a horocycle
                        EA[0]=EA[1];
                        EA[1]=t;
                    }
                    *cw=((dat)>>3)&511;
                    int di=int(paw[*cw])+sizeOfPow*(dat>>12)&1023;
                    si = di;
                    sp1[0]=paw[di];sp1[1]=paw[di+1];sp1[2]=paw[di+2];
                    sp2[0]=paw[di+3];sp2[1]=paw[di+4];sp2[2]=paw[di+5];
                    d+=ped-cid;
                }
            } else {//to do: remove mirror from here
                pos[0]=d;pos[1]=0;pos[2]=0;
                pos[0]-=off[0];pos[1]-=off[1];
                rotate(pos[0],pos[1],rot[0],-rot[1]);//reversing rotations, which is easy in E2
                pos2r[1]*=mirror;
                pos2[0]=pos2r[0]*pos2r[2]+d;pos2[1]=pos2r[1]*pos2r[2];
                pos2[0]-=off[0];pos2[1]-=off[1];pos2[2]=0;
                rotate(pos2[0],pos2[1],rot[0],-rot[1]);
                //printf("%f %f\n",pos2[0],pos2[1]);
                ref[0]=dx+d;ref[1]=-dy*mirror;ref[2]=0;
                ref[0]-=off[0];ref[1]-=off[1];
                rotate(ref[0],ref[1],rot[0],-rot[1]);
                d=0;
                //if(iterations>9) bounce=1;

                //TO DO: Make wall checks its own function?
                if(!bounce) {
                    //*
                    //
                    double RL[3], pos2r[3], posc[2];
                    double closestpt[2];
                    int i = (int)paw[*cw];
                    bool inwall = 0, inportal = 0;
                    while(i < paw[*cw + 1] && !inwall){
                        bool valid = 0;
                        double cp[2];
                        double p1[] = {paw[i], paw[i + 1]}, p2[] = {paw[i + 3],paw[i + 4]};
                        p2[0] -= p1[0]; p2[1] -= p1[1];
                        posc[0] = pos[0] - p1[0]; posc[1] = pos[1] - p1[1];
                        cp[0] = posc[0]; cp[1] = posc[1];
                        double rl2 = p2[0]*p2[0] + p2[1]*p2[1];//radius/length^2
                        double rc;
                        type = ((int)paw[i + 7])&7;
                        if( (type&1) == 0){
                            rc = sqrt(rl2/(posc[0]*posc[0] + posc[1]*posc[1]) );//r conversion factor
                            cp[0]*=rc;cp[1]*=rc;
                            double al = paw[i + 6], ala = rl2*al;//angle limit (2-2cos(a)) and adjusted angle limit
                            if(al==4 || ( (cp[0] - p2[0])*(cp[0] - p2[0]) + (cp[1] - p2[1])*(cp[1] - p2[1]) < ala) ){
                                valid = 1;
                                //...
                            } else if(type > 3){
                                double c = 1 - al*0.5, s = sqrt(1 - c*c);
                                if(posc[1]*p2[0] - posc[0]*p2[1] < 0) s = -s;
                                cp[0] = p2[0]; cp[1] = p2[1];
                                rotate(cp[0],cp[1],c,s);
                            }
                        } else {
                            t=posc[0]*p2[0] + posc[1]*p2[1];
                            valid = 1;
                            if(t < 0) t = 0;
                            if(t > rl2) t = rl2;
                            cp[0] = p2[0]*t/rl2;
                            cp[1] = p2[1]*t/rl2;
                        }
                        t = (posc[0] - cp[0])*(posc[0] - cp[0]) + (posc[1] - cp[1])*(posc[1] - cp[1]);
                        inwall = t < playerWallCollisionr*playerWallCollisionr && type > 3;
                        if( inwall || (t < pr*pr && !inportal && valid) ){
                            RL[2] = sqrt(t);
                            if(!inwall){
                                inportal = 1;
                                if(t < ped*ped) inwall = 1;
                                else {
                                    if(type == 0){
                                        RL[0] = cp[0]; RL[1] = cp[1];
                                        rotate(RL[0],RL[1],p2[0],-p2[1]);
                                        RL[0] /= rl2; RL[1] /=rl2;
                                        if(rc > 1) RL[2] = -RL[2];
                                    } else {
                                        rl2 = sqrt(rl2);
                                        RL[1] = rl2/2 - (cp[0]*p2[0] + cp[1]*p2[1])/rl2;
                                        if(posc[1]*p2[0] - posc[0]*p2[1] < 0) RL[2] = -RL[2];
                                    }
                                    cp[0] -= posc[0];cp[1] -= posc[1];
                                    double pos2c[2] = {pos2[0] - pos[0], pos2[1] - pos[1]};
                                    double pos2r[3] = {-(pos2c[0]*cp[0] + pos2c[1]*cp[1])/pr/RL[2], (pos2c[0]*cp[1] - pos2c[1]*cp[0])/pr/RL[2], pr};
                                    double pos1[3], pos3[3];//this is jank
                                    createDuplicate(RL,RL[2],pos2r,i,&duppw,pos1,pos3);
                                    copypt(pos1,duppl);
                                    copypt(pos3,duppl2);
                                    iopip=i;//make this entity specific later
                                }
                            }
                            if(inwall){
                                d = sqrt(t);
                                cp[0] -= posc[0];cp[1] -= posc[1];
                                double refc[2] = {ref[0] - pos[0], ref[1] - pos[1]};
                                dx = -(refc[0]*cp[0] + refc[1]*cp[1])/d;
                                dy = -(refc[0]*cp[1] - refc[1]*cp[0])/d;
                                if(inportal) d = 1.01*ped - d;
                                else d = 1.01*playerWallCollisionr - d;
                                mirror=1;
                                type=0;si=-1;
                                sp1[0]=pos[0];
                                sp1[1]=pos[1];
                                sp2[0]=ref[0];
                                sp2[1]=ref[1];
                                EL[0]=dx;EL[1]=-dy;
                            }
                            //...
                        }
                        //
                        i+=sizeOfPow;
                    }
                }//end of collision checks
                //else d=0;
;
                //if(i==paw[*cw+1]) d=0;
            }//end of if there are no portals




            //[-x  y][-1]   [x]
            //[-y -x][ 0] = [y]

        } else if(iaunit > 0){//s2
            double dc=cos(d),ds=sin(d);
            double rot[3][3];
            double rc,rs,r;
            if(type < 4){
                if(si > 0){
                    rc = dot(sp1,sp2);
                    rs = sqrt(1 - rc*rc);
                    s2matto(si,rot);
                } else {
                    s2matto(sp1,rot);
                    matxpt(rot,sp2);
                    r=sqrt(sp2[0]*sp2[0]+sp2[1]*sp2[1]);
                    sp2[0]/=r;sp2[1]/=-r;
                    //way more effecient than matrix x matrix multiplication, which will be avoided at all costs on the CPU side.
                    rotXY(rot,sp2[0],sp2[1]);
                    //end of axis aligning sp2, well making the transformation that would.
                }
                rotXY(rot,EL[0],EL[1]);
                if(si > 0){
                    rotXZ(rot,rc,rs);
                    rotXY(rot,EA[0],EA[1]);
                } else {
                    copypt(pos2,pos2r);
                    matxpt(rot,pos2r);
                    r=sqrt(pos2r[0]*pos2r[0]+pos2r[1]*pos2r[1]);
                    pos2r[0]/=r;pos2r[1]/=r;
                    pos2r[2]=pr;
                }
            } else {
                copymat(rotmstart,rot);
            }
            //checks for portals
            int i=(int)paw[*cw],cii=-1;
            double cip1[3],cip2[3],cip[3];
            double cidr=s2disrank(dc,ds);
            while(i<(int)paw[*cw+1]){
                double p1[3]={paw[i],paw[i+1],paw[i+2]}, p2[3]={paw[i+3],paw[i+4],paw[i+5]};
                rc=dot(p1,p2);
                rs=1-rc*rc;//actually rs^2, just reusing space
                matxpt(rot,p1);
                matxpt(rot,p2);
                double y1y1=p1[1]*p1[1];
                if(y1y1 < rs){
                    double cipp[3], cipm[3];//closest intersection point plus and minus
                    double s=sqrt(rs-y1y1), D=1-y1y1;//angle limit, (cos(a)-1), decreases with angle
                    double al=paw[i+6], ala=1+rs*al, disrank;
                    cipp[0]=(p1[0]*rc+p1[2]*s)/D;
                    cipm[0]=(p1[0]*rc-p1[2]*s)/D;
                    cipp[2]=(p1[2]*rc-p1[0]*s)/D;
                    cipm[2]=(p1[2]*rc+p1[0]*s)/D;
                    cipp[1]=0;cipm[1]=0;

                    disrank=s2disrank(cipp[2],cipp[0]);
                    if( ( al==-2 || dot(cipp,p2)>ala ) && disrank<cidr && (si!=i || cipp[2] < cipm[2]) ){
                        cii=i;
                        copypt(p1,cip1);
                        copypt(p2,cip2);
                        copypt(cipp,cip);
                        cidr=disrank;
                    }
                    disrank=s2disrank(cipm[2],cipm[0]);
                    if( ( al==-2 || dot(cipm,p2)>ala ) && disrank<cidr && (si!=i || cipm[2] < cipp[2]) ){
                        cii=i;
                        copypt(p1,cip1);
                        copypt(p2,cip2);
                        copypt(cipm,cip);
                        cidr=disrank;
                    }

                }
                i+=sizeOfPow;
            }//end of portal checks*/
            //
            if(cii > -1){
                int dat=(int)paw[cii+7];
                type = dat&7;
                if(type > 3){
                    si = cii;
                    d -= arctan(cip[0],cip[2]);
                    rotate(cip1[0],cip1[2],cip[2],cip[0]);
                    rotXZ(rot,cip[2],cip[0]);
                    double r2 = cip1[0]*cip1[0] + cip1[1]*cip1[1];
                    double c = 2*cip1[1]*cip1[1]/r2 - 1, s = -2*cip1[0]*cip1[1]/r2;
                    rotXY(rot,c,-s);
                    rotate(pos2r[0],pos2r[1],c,-s);
                    copymat(rot,rotmstart);
                    //...
                } else {
                    char side;
                    if(((dat>>22)&1)==1) side=-1;
                    else side=1;
                    if(((dat>>23)&1)==1) {mirror=-1;props[0]=-props[0];}
                    else mirror=1;
                    *cw=(dat>>3)&511;
                    int di=(int)(paw[*cw])+sizeOfPow*(dat>>12)&1023;
                    vec3(sp1, paw[di],paw[di+1],paw[di+2]);
                    vec3(sp2, paw[di+3],paw[di+4],paw[di+5]);
                    si=di;
                    d -= arctan(cip[0],cip[2]);
                    {
                    double cosr2 = dot(cip1,cip2); cosr2 *= cosr2;
                    double cosE = (dot(cip2,cip) - cosr2)/(1 - cosr2);
                    double sinE = -safe_sqrt(1 - cosE*cosE)*mirror*side;
                    //getting the sign of sin, dot(cross(p2,p1),p)
                    //to do: make this its own function later
                    double cross[3];
                    cross[0] = cip1[1]*cip2[2] - cip1[2]*cip2[1];
                    cross[1] = -cip1[0]*cip2[2] + cip1[2]*cip2[0];
                    cross[2] = cip1[0]*cip2[1] - cip1[1]*cip2[0];
                    if(dot(cip,cross) < 0) sinE = -sinE;
                    EL[0] = cosE;
                    EL[1] = sinE;
                    }
                    double tm[3][3];//temporary matrix
                    rotate(cip1[0],cip1[2],cip[2],cip[0]);
                    r=sqrt(cip1[0]*cip1[0]+cip1[1]*cip1[1]);
                    cip1[0]/=r;cip1[1]/=r;
                    EA[0]=-cip1[0]*side;
                    EA[1]=-cip1[1]*mirror*side;
                }
            } else {
                double roti[3][3];
                transpose(rot,roti);
                vec3(pos, ds,0,dc);
                s2dadtopt(pos2r[2],pos2r[0],pos2r[1]*mirror,pos2);
                vec3(ref, dx,-dy*mirror,0);
                rotate(pos2[2],pos2[0],dc,ds);
                rotate(ref[2],ref[0],dc,ds);
                matxpt(roti,pos);
                matxpt(roti,pos2);
                matxpt(roti,ref);
                backOnSphere(pos);//stopping errors from accumilating beyond a point
                d=0;
                if(!bounce){
                    double cp[3], closestpt[3], ip1[3], ip2[3];
                    double cosODis;
                    i = (int)paw[*cw];
                    int indexOfIntersection = -1;
                    char inwall = 0, inportal = 0;
                    while(i<(int)paw[*cw+1] && inwall==0){
                        bool valid = 1;
                        char checkifinportal = 1;
                        double p1[3]={paw[i],paw[i+1],paw[i+2]}, p2[3]={paw[i+3],paw[i+4],paw[i+5]};
                        double rc = dot(p1,p2), //cos(r)
                        rss = 1 - rc*rc,        //sin(r)^2
                        r2c = dot(pos,p1),      //cos(r from entity to center of circle)
                        A = safe_sqrt( rss/(1 - r2c*r2c)),
                        B = rc - A * r2c;
                        cp[0] = A*pos[0] + B*p1[0];
                        cp[1] = A*pos[1] + B*p1[1];
                        cp[2] = A*pos[2] + B*p1[2];
                        char iswall = ( (int)paw[i + 7])&4;
                        double al = paw[i + 6], ala = 1 + rss*al;
                        if( (al > -2.0 && dot(cp,p2) < ala || r2c > cosplayerWallCollisionr)  ){
                            if(iswall > 0){
                                double endpt1[3] = {paw[i + 8],paw[i + 9],paw[i + 10]},
                                endpt2[3] = {paw[i + 11],paw[i + 12],paw[i + 13]};
                                double dot1 = dot(pos,endpt1), dot2 = dot(pos,endpt2);
                                if(dot1>dot2) copypt(endpt1,cp);
                                else copypt(endpt2,cp);
                            } else valid = 0;
                        }
                        double cosDisToCP = dot(pos,cp);
                        if( (iswall > 0 && cosDisToCP > cosplayerWallCollisionr) || (iswall == 0 && valid && cosDisToCP > prc && inportal == 0) ){
                            cosODis = cosDisToCP;
                            if(iswall > 0) {
                                inwall = 1;
                                inportal = 0;
                            } else {
                                inportal = 1;
                                copypt(p1,ip1);
                                copypt(p2,ip2);
                            }
                            copypt(cp,closestpt);
                            indexOfIntersection = i;
                        }
                        //...
                        i += sizeOfPow;
                    }
                    if(cosODis > 1) cosODis = 1;
                    if(cosODis < -1) cosODis = -1;
                    if(cosODis > cos(ped)) inwall = 1;
                    if(indexOfIntersection > 0){
                        if(inportal > 0 && !inwall){
                            double posc[3], pos2c[3], relLoc[3];
                            copypt(pos,posc);
                            copypt(pos2,pos2c);
                            double BC[3], AB[3];
                            cross(pos2c,posc,AB);
                            AB[0] /= prs; AB[1] /= prs; AB[2] /= prs;
                            cross(posc,closestpt,BC);
                            backOnSphere(BC);
                            double cr2 = dot(ip1,ip2);
                            if(dot(ip1,posc) > cr2) relLoc[2] = -1;
                            else relLoc[2] = 1;
                            double facingCos = dot(AB,BC)*relLoc[2];
                            double facingSin = -safe_sqrt(1 - facingCos*facingCos)*sign(dot(BC,pos2c))*relLoc[2];
                            double temppos2r[] = {facingCos,facingSin,pr};
                            cr2 *= cr2;
                            relLoc[0] = (dot(closestpt,ip2) - cr2)/(1 - cr2);
                            relLoc[1] = safe_sqrt(1 - relLoc[0]*relLoc[0])*s2side(ip1,ip2,closestpt);
                            relLoc[2] *= acos(cosODis);
                            double pos1[3], pos3[3];//temp sotrage due to type conversion issues
                            createDuplicate(relLoc,relLoc[2],temppos2r,indexOfIntersection,&duppw,pos1,pos3);
                            copypt(pos1,duppl);
                            copypt(pos3,duppl2);
                            iopip = indexOfIntersection;//make this entity specific later
                        }
                        if(inwall > 0){
                            if(inportal > 0) d = 1.01*ped - acos(cosODis);
                            else d = 1.01*playerWallCollisionr - acos(cosODis);
                            double refc[3];
                            copypt(ref,refc);
                            s2matto(pos,rot);
                            matxpt(rot,refc);
                            matxpt(rot,closestpt);
                            rotate(closestpt[0],closestpt[1],refc[0],-refc[1]);
                            double t = sqrt(closestpt[0]*closestpt[0] + closestpt[1]*closestpt[1]);
                            dx = -closestpt[0]/t;
                            dy = -closestpt[1]/t;
                            EL[0] = dx;
                            EL[1] = -dy;
                            si = -1;
                            type = 0;
                            mirror = 1;
                            copypt(pos,sp1);
                            copypt(ref,sp2);
                        }
                    }
                }//*/
            }
        }//end of s2
        else if(iaunit < 0){//h2
            double dc = exp(d), ds = 1/dc;
            dc = (dc + ds)/2; ds = dc - ds;
            //^cosh and sinh of d, done weird to remove redundant calculations
            double k;
            double rot[3][3];
            if(type < 4){
                k = lidot(sp1,sp2);
                if(type > 0){
                    double tv[3];
                    copypt(sp1,tv);
                    copypt(sp2,sp1);
                    copypt(tv,sp2);
                }
                h2matto(sp1,rot);
                matxpt(rot,sp2);
                double r = sqrt(sp2[0]*sp2[0] + sp2[1]*sp2[1]);
                sp2[0] /= r; sp2[1] /= r;
                rotXY(rot, sp2[0],-sp2[1]);
                switch(type){
                    case 0:
                        rotXY(rot, EL[0],EL[1]);
                        if(si < 0){
                            copypt(pos2,pos2r);
                            matxpt(rot,pos2r);
                            r = sqrt(pos2r[0]*pos2r[0] + pos2r[1]*pos2r[1]);
                            pos2r[0] /= r;pos2r[1] /= r;
                            pos2r[2] = pr;
                        } else lorenzXZ(rot,k,-r);
                        break;
                    case 1:
                        idealrot(rot,EL[1]);
                        break;
                    case 2:
                        lorenzXZ(rot,r,-k);
                        lorenzYZ(rot,EL[0],-EL[1]);
                        lorenzXZ(rot,r,k);
                        break;
                    default:
                        printf("invalid type\n");
                }
                if(si > 0) rotXY(rot,EA[0],EA[1]);
            } else {
                copymat(rotmstart,rot);
            }
            double cip1[3], cip2[3], cip[3];
            cip[0] = ds;
            int cii = -1;
            int i = int(paw[*cw]);
            while(i < (int)paw[*cw + 1]){
                double p1[3] = {paw[i],paw[i + 1],paw[i + 2]},
                p2[3] = {paw[i + 3],paw[i + 4],paw[i + 5]};
                double limit = paw[i + 6];
                type = int(paw[i + 7])&3;
                k = lidot(p1,p2);
                matxpt(rot,p1);
                matxpt(rot,p2);
                double o = sqrt(k*k - p1[2]*p1[2] + p1[0]*p1[0]);
                double tR = (k + o)/(p1[2] - p1[0]), itR = 1/tR;
                double tL = (k + o)/(p1[2] + p1[0]), itL = 1/tL;
                double cipR[3], cipL[3];
                cipR[1] = cipL[1] = 0;
                cipR[0] = (tR - itR)/2; cipR[2] = cipR[0] + itR;
                cipL[0] = (itL - tL)/2; cipL[2] = cipL[0] + tL;
                bool risvalid = tR > 0 && limit > lidot(p2,cipR);
                bool lisvalid = tL > 0 && limit > lidot(p2,cipL);
                if(risvalid && 0 < cipR[0] && cipR[0] < cip[0] && (i != si || (lisvalid && cipL[2] < cipR[2]))){
                    copypt(cipR,cip);
                    cii = i;
                }
                if(lisvalid && 0 < cipL[0] && cipL[0] < cip[0] && (i != si || (risvalid && cipR[2] < cipL[2]))){
                    copypt(cipL,cip);
                    cii = i;
                }
                if(i == cii){
                    copypt(p1,cip1);
                    copypt(p2,cip2);
                }
                i += sizeOfPow;
            }
            if(cii > 0){
                int dat = (int)paw[cii + 7];
                type = dat&7;
                if(type > 3){
                    si = cii;
                    d -= acosh(cip[2]);
                    lorenz(cip1[0],cip1[2],cip[2],-cip[0]);
                    lorenzXZ(rot,cip[2],-cip[0]);
                    double r2 = cip1[0]*cip1[0] + cip1[1]*cip1[1];
                    double c = 2*cip1[1]*cip1[1]/r2 - 1, s = -2*cip1[0]*cip1[1]/r2;
                    rotXY(rot,c,-s);
                    rotate(pos2r[0],pos2r[1],c,-s);
                    copymat(rot,rotmstart);
                } else {
                    d -= acosh(cip[2]);
                    double side = 1;
                    if(((dat>>22)&1)==1) side = -1;
                    mirror = 1;
                    if(((dat>>23)&1)==1) {mirror=-1;props[0]=-props[0];}//
                    *cw = (dat>>3)&511;
                    int di=(int)(paw[*cw])+sizeOfPow*(dat>>12)&1023;
                    vec3(sp1, paw[di],paw[di+1],paw[di+2]);
                    vec3(sp2, paw[di+3],paw[di+4],paw[di+5]);
                    si = di;
                    k = lidot(cip1,cip2);
                    double k2 = k*k;
                    double cL = lidot(cip,cip2), sL = 1;
                    switch(type){
                        case 0:
                            cL = (cL - k2)/(1 - k2);
                            sL = 1 - cL*cL;
                            break;
                        case 1:
                            sL = 2*cL - 2;
                            break;
                        case 2:
                            cL = (cL + k2)/(k2 + 1);
                            sL = cL*cL - 1;
                            break;
                        default:
                            printf("invalid type\n");
                    }
                    sL = safe_sqrt(sL)*h2side(cip1,cip2,cip);
                    EL[0] = cL;
                    EL[1] = sL*side*mirror;
                    lorenz(cip1[0],cip1[2],cip[2],-cip[0]);
                    if(type > 0){
                        cip1[0] = -cip1[0];
                        cip1[1] = -cip1[1];
                    }
                    double r = sqrt(cip1[0]*cip1[0] + cip1[1]*cip1[1]);
                    EA[0] = -cip1[0]/r*side;
                    EA[1] = -cip1[1]/r*side*mirror;
                }
            } else {
                //...
                //bug: camera not placed right when coming in from portal
                double roti[3][3];
                h2invert(rot,roti);
                vec3(pos, ds,0,dc);
                h2dadtopt(pos2r[2],pos2r[0],pos2r[1]*mirror,pos2);
                vec3(ref, dx,-dy*mirror,sqrt(2));
                lorenz(ref[0],ref[2],dc,ds);
                lorenz(pos2[0],pos2[2],dc,ds);
                matxpt(roti,pos);
                matxpt(roti,pos2);
                matxpt(roti,ref);
                backOnHyperboloid(pos);
                backOnHyperboloid(pos2);
                backOnHyperboloid(ref);
                d = 0;
                if(!bounce){
                    int i = paw[*cw];
                    double posc[3];
                    copypt(pos,posc);
                    bool inwall = false;
                    bool inportal = false;
                    double coshdis;
                    double closestpoint[3];
                    while(i < paw[*cw + 1] && !inwall){
                        bool valid = 0;
                        double p1[3] = {paw[i],paw[i + 1],paw[i + 2]},
                        p2[3] = {paw[i + 3],paw[i + 4],paw[i + 5]};
                        double limit = paw[i + 6];
                        type = int(paw[i + 7])&3;
                        bool iswall = (int(paw[i + 7])&4) > 0;
                        k = lidot(p1,p2);
                        double id1 = lidot(posc,p1);
                        double id2 = lidot(p1,p1);
                        double A = id1*id1 - id2;
                        double B = sqrt( (k*k - id2)/A );
                        double C;
                        if(type == 1) C = (id1*id1 - k*k)/(2*A*k);
                        else C = (k - B*id1)/id2;
                        double cp[3] = {B*posc[0] + C*p1[0],B*posc[1] + C*p1[1],B*posc[2] + C*p1[2]};
                        if(lidot(cp,p2) < limit){
                            valid = 1;
                            //...
                        } else {
                            int off;
                            if(h2side(p1,p2,cp) > 0) off = 8;
                            else off = 11;
                            vec3(cp, paw[i + off],paw[i + off + 1],paw[i + off + 2]);
                        }
                        double c = lidot(cp,posc);
                        if((c < prch && !inportal && valid) || (iswall && c < coshplayerWallCollisionr)){
                            copypt(cp, closestpoint);
                            coshdis = c;
                            if(coshdis < 1) coshdis = 1;
                            if(iswall){
                                inwall = true;
                                inportal = false;
                            } else {
                                inportal = true;
                                if(coshdis > cosh(ped)){
                                        double relLoc[3];
                                        double k2 = k*k;
                                        double cL = lidot(cp,p2), sL = 1;
                                        switch(type){
                                            case 0:
                                                cL = (cL - k2)/(1 - k2);
                                                sL = 1 - cL*cL;
                                                break;
                                            case 1:
                                                sL = 2*cL - 2;
                                                break;
                                            case 2:
                                                cL = (cL + k2)/(k2 + 1);
                                                sL = cL*cL - 1;
                                                break;
                                            default:
                                                printf("invalid type\n");
                                        }
                                        sL = -safe_sqrt(sL)*h2side(p1,p2,posc);
                                        double pos2c[3], AB[3], BC[3];
                                        copypt(pos2,pos2c);
                                        licross(pos2c,posc,AB);
                                        AB[0] /= prsh; AB[1] /= prsh; AB[2] /= prsh;
                                        licross(posc,cp,BC);
                                        double t = sqrt(-lidot(BC,BC));
                                        BC[0] /= t; BC[1] /= t; BC[2] /= t;
                                        if(lidot(posc,p1) > k) relLoc[2] = 1;//checking which side of the object I am on for format reasons
                                        else relLoc[2] = -1;
                                        double facingCos = -lidot(AB,BC)*relLoc[2];
                                        double facingSin = safe_sqrt(1 - facingCos*facingCos)*sign(lidot(BC,pos2c))*relLoc[2];
                                        double temppos2r[] = {facingCos,facingSin,pr};
                                        double cr2 = coshdis*coshdis;
                                        relLoc[0] = cL;
                                        relLoc[1] = sL;
                                        relLoc[2] *= acosh(coshdis);
                                        double pos1[3], pos3[3];//temp sotrage due to type conversion issues
                                        createDuplicate(relLoc,relLoc[2],temppos2r,i,&duppw,pos1,pos3);
                                        copypt(pos1,duppl);
                                        copypt(pos3,duppl2);
                                        iopip = i;//TODO make this entity specific later
                                } else {
                                    inwall = true;
                                }
                            }
                        }
                        i += sizeOfPow;
                    }
                    if(inwall){
                        double rot[3][3];
                        h2matto(posc,rot);
                        double refc[3];
                        copypt(ref,refc);
                        matxpt(rot,refc);
                        matxpt(rot,closestpoint);
                        double r = sqrt(refc[0]*refc[0] + refc[1]*refc[1]);
                        rotate(closestpoint[0],closestpoint[1],refc[0]/r,-refc[1]/r);
                        //printf("%.99lf\n",coshdis);
                        if(coshdis < 1) coshdis = 1;
                        r = sqrt(closestpoint[0]*closestpoint[0] + closestpoint[1]*closestpoint[1]);
                        if(inportal) d = 1.01*ped - acosh(coshdis);
                        else d = 1.01*playerWallCollisionr - acosh(coshdis);
                        //printf("%.99lf\n",d);
                        //printf("%.99lf\n",r);
                        dx = -closestpoint[0]/r;
                        dy = -closestpoint[1]/r;
                        //printf("%lf\t%lf\t%lf\n",dx,dy,d);
                        EL[0] = dx;
                        EL[1] = -dy;
                        si = -1;
                        type = 0;
                        mirror = 1;
                        copypt(pos,sp1);
                        copypt(ref,sp2);
                        //d = 0;
                    }
                }//end wall and dup stuff
            }
        }//end of h2*/
    }
    if(ogcw == *cw){
        double posc[3];
        copypt(pos,posc);
        posc[0] -= ogpos[0];posc[1] -= ogpos[1];posc[2] -= ogpos[2];
        if(posc[0]*posc[0] + posc[1]*posc[1] + posc[2]*posc[2] < ped*ped){//reverting if the move was less than a rounding error. This is to fix a bug with wall corners in S2
            copypt(ogpos,pos);
            copypt(ogpos2,pos2);
            copypt(ogref,ref);
        }
    }
    //std::cout<<iterations<<"\n";
    if(iterations >= maxits) std::cout<<"failsafe trigger\n";
    //if(*cw==0 &&  pos[0]*pos[0] + (pos[1]-2)*(pos[1]-2) > 1.25*1.25 ) printf("error?\n");
}


/*class e2chunk{
public:
    int x;
    int y;
    std::vector<std::vector<double>> geo;
};

class h2chunk{
public:
    std::vector<char> address;
    std::vector<std::vector<double>> geo;
};*/


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

int main(){
    //taking in account that the right stick is mapped differently on Windows than on Linux
    int rsxi, rsyi;

    #ifdef __GNUC__
        rsxi = 3;
        rsyi = 4;
    #endif
    #ifdef _WIN64
        rsxi = 2;
        rsyi = 3;
    #endif
    /*{//psudosphere case b
        double r = 2;
        double sr = sinh(r), cr = cosh(r), x2 = sinh(r/2), z2 = cosh(r/2);
        double a = 2*pi/cr, c1 = cos(a), s1 = sin(a), c2 = cos(a/2), s2 = sin(a/2);
        double limit = cr*cr - sr*sr*c2;
        std::vector<std::vector<float>> world0 = {
            {0,1,0, x2,0,z2, z2,mkdest(2,0,1,0,0)},
            {-s1,c1,0, x2*c1,x2*s1,z2, z2,mkdest(2,0,0,0,0)},
            {0,0,1, sr*c2,sr*s2,cr, limit,mkdest(0,1,2,1,0)}
        };
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world1 = {
            {0,1,0, x2,0,z2, z2,mkdest(2,1,1,0,0)},
            {-s1,c1,0, x2*c1,x2*s1,z2, z2,mkdest(2,1,0,0,0)},
            {0,0,1, sr*c2,sr*s2,cr, limit,mkdest(0,0,2,1,0)}
        };
        pawbuffer.push_back(world1);
        worldCurvatures.push_back(-1);
        pw = 0; plp[0] = 1; duppw = -1;
        vec3(pl, 1,1,sqrt(3)); vec3(camRef, 1,0,sqrt(2));
        vec3(pl2, prsh,0,prch);
    }//*/
    /*{//psudosphere case c
        double d = 1.5, ds = sinh(d), dc = cosh(d), l = pi/ds, sl = sinh(l), lc = cosh(l);
        double limit = dc*dc*lc - ds*ds;
        std::vector<std::vector<float>> world0={
        {lc,0,sl, sl,0,lc, dc,mkdest(2,0,2,0,0)},
        {0,1.00,0, 0,ds,dc, limit,6+0*mkdest(2,0,3,0,0)},
        {lc,0,-sl, -sl,0,lc, dc,mkdest(2,0,0,0,0)},
        {0,1.00,0, 0,-ds,dc, limit,6+0*mkdest(2,0,1,0,0)}
        };
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(-1);
        pw = 0; plp[0] = 1; duppw = -1;
        vec3(pl, 0,0,1); vec3(camRef, 1,0,sqrt(2));
        vec3(pl2, prsh,0,prch);
    }//*/

    /*{//S2 portal cutoff test
        std::vector<std::vector<float>> world0={
        {1,0,0, isr2,0,isr2, -2,4+0*mkdest(0,0,1,0,0)},
        {1,0,0, -isr2,0,isr2, -2,4+0*mkdest(0,0,0,0,0)}
        };//d*d/2 + 1 = -L for horocycles
        //std::vector<std::vector<float>> world0={{0,3,sqrt(8), 0,1,sqrt(2), 3,6}};
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(1);
        pw = 0; plp[0] = 1; duppw = -1;//signals that the player isn't in a portal
        vec3(pl, 0,0,1); vec3(camRef, 1,0,0);
        //
        //vec3(pl, 0,0,0);vec3(pl2, pl[0],pl[1]+pr,0);vec3(camRef, pl[0]+1,pl[1],0);pw=1;
    }//*/

    /*{//H2 portal cutoff test
        std::vector<std::vector<float>> world0={
        {2,0,sqrt(4), 1,0,sqrt(2), 2,mkdest(1,0,1,0,1)},
        {-2,0,sqrt(4), -1,0,sqrt(2), 2,mkdest(1,0,0,0,1)}
        };//d*d/2 + 1 = -L for horocycles
        //std::vector<std::vector<float>> world0={{0,3,sqrt(8), 0,1,sqrt(2), 3,6}};
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(-1);
        pw = 0; plp[0] = 1; duppw = -1;//signals that the player isn't in a portal
        vec3(pl, 0,0,1); vec3(camRef, 1,0,sqrt(2));
        vec3(pl2, prsh,0,prch);
        //
        //vec3(pl, 0,0,0);vec3(pl2, pl[0],pl[1]+pr,0);vec3(camRef, pl[0]+1,pl[1],0);pw=1;
    }//*/

    /*{//h2 e2 test
        std::vector<std::vector<float>> world0={
            {0,2,sqrt(4), 0,0,1, 3,mkdest(1,1,0,0,0)},//d*d/2 + 1 = -L for horocycles
            {0,-3,sqrt(8), 0,-1,sqrt(2), 3,6}
        };
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world1={{2,3,0, -2,3,0, 4,mkdest(1,0,0,0,0)}};
        pawbuffer.push_back(world1);
        worldCurvatures.push_back(0);
        pw = 0; plp[0] = 1; duppw = -1;//signals that the player isn't in a portal
        vec3(pl, 0,0,1); vec3(camRef, 1,0,sqrt(2));
        vec3(pl2, prsh,0,prch);
        //
        //vec3(pl, 0,0,0);vec3(pl2, pl[0],pl[1]+pr,0);vec3(camRef, pl[0]+1,pl[1],0);pw=1;
    }//*/
    //NOTE: isr2 = 1/sqrt(2), sr2 = sqrt(2), I just got tired of typing it in testing

    {//test stage
        std::vector<std::vector<float>> world0={
            {0,2,0, 0.5,2,0, 2,mkdest(0,1,0,1,0)},//portal, type 0 (arc), linked to world1, index 0, other side = 1, mirror = 0
            {2,4,0, 0,4,0, 0,mkdest(1,0,1,0,1)},//portal, type 1(E2 line), linked to world 0, index 1 (itself), side = 0, mirror = 1
            {0,-5,0, 0,-4,0, 2,mkdest(0,0,3,0,0)},
            {3,-2,0, 2,-2,0, 2,mkdest(0,0,2,0,0)},
            {-2,5,0, -2,1,0, 0,mkdest(1,2,0,0,0)},
            {1,1,0, 2,1,0, 0,5},//E2 line wall
            {2,2,0, 3,2,0, 2,4},//arc wall
        };
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(0);//euclidian
        /////////////////////////
        std::vector<std::vector<float>> world1={
        {0,0,1, 0,0.5,0.866025403784, -1,mkdest(0,0,0,1,0)}, //portal, type 0 (arc), linked to world0, index 0, other side = 1, mirror = 0
        {0,0,1, 0,1,0, -0.5,4}, //arc wall
        {-1,0,0, 0,isr2,-isr2, -0.292893218813,4}//arc wall
        };
        pawbuffer.push_back(world1);
        worldCurvatures.push_back(1);//spherical
        /////////////////////
        //coments WIP
        std::vector<std::vector<float>> world2={
            {0,2,sqrt(4), 0,0,1, 3,mkdest(1,0,4,0,0)},//d*d/2 + 1 = -L for horocycles
            {0,-3,sqrt(8), 0,-1,sqrt(2), 3,6}//hyperciclic wall
        };
        pawbuffer.push_back(world2);
        worldCurvatures.push_back(-1);//hyperbolic
        //
        pw = 0;//player in world 0
        vec3(pl, 0,0,0);//player at 0,0,0
        vec3(pl2, pl[0],pl[1]+pr,0);//placing the default facing of the player
        vec3(camRef, pl[0]+1,pl[1],0);
        plp[0] = 1;//setting this to negative 1 mirrors
        duppw = -1;//signals that the player isn't in a portal
    }//*/



    /*{
        //std::vector<std::vector<float>> world={ {0,1,0, 0,isr2,isr2, -0.5, mkdest(0,0,1,0,0)} , {0,0,1, -isr2,0,isr2, -0.5, mkdest(0,0,0,0,0)}};//s2
        //std::vector<std::vector<float>> world={ {0,0,1, isr2,0,isr2, -0.5, mkdest(0,0,1,1,1)} , {0,0,1, -isr2,0,isr2, -0.5, mkdest(0,0,0,1,1)}, {0,0,1, 0,1,0, -1.05,4} };//s2
        //std::vector<std::vector<float>> world={{0,0,1, 0,isr2,isr2, -0.25,4},{0,0,1, -isr2,0,isr2, -0.25,4}};

        //std::vector<std::vector<float>> world={{0,0,0,3,0,0,2,mkdest(0,0,1,1,1)},{0,0,0,-3,0,0,2,mkdest(0,0,0,1,1)}};

        //std::vector<std::vector<float>> world={{4,0,0,2,0,0,1,mkdest(0,0,1,0,0)},{-1*isr2,1*isr2,0,-3*isr2,3*isr2,0,1,mkdest(0,0,0,0,0)},{-2,-3,0,-4,-1,0,0,mkdest(1,1,0,0,0)},{-1,-1,0,-1,-2,0,4,mkdest(0,1,1,1,0)},{1,-4,0,1,-2,0,3,4}};

        //type, world,index, side?, mirror?
        //std::vector<std::vector<float>> world={{-1*isr2,1*isr2,0,-3*isr2,3*isr2,0,2,mkdest(0,0,0,0,1)}};

        //std::vector<std::vector<float>> world={{2,0,0,4,0,0,0,mkdest(1,0,1,0,1)},{-1*isr2,1*isr2,0,-3*isr2,3*isr2,0,0,mkdest(1,0,0,0,1)},{-1,-3,0,1,-3,0,3,5},{2,2,0,1,1,0,2,mkdest(0,0,4,1,0)},{-3,-2,0,-2,-1,0,2,mkdest(0,0,3,1,0)},{1,-1.5,0,2.5,-1.5,0,2.01,4}};

        //std::vector<std::vector<float>> world={{2,0,0,4,0,0,0,mkdest(1,0,1,0,1)},{-2,0,0,-2,2,0,0,mkdest(1,0,0,0,1)},{-1,-3,0,1,-3,0,3,5},{1,1,0,2,2,0,2,mkdest(0,0,4,1,0)},{-3,-2,0,-2,-1,0,2,mkdest(0,0,3,1,0)},{1,-1.5,0,2.5,-1.5,0,2.01,4},{3,3,0,4,4,0,2,mkdest(0,0,7,0,0)},{-4,-3,0,-3,-2,0,2,mkdest(0,0,6,0,0)}};

        //std::vector<std::vector<float>> world={{1,-1.5,0,2.5,-1.5,0,2.01,5}};

        std::vector<std::vector<float>> world={
        {-1,1.5,0,0,1.5,0,0,mkdest(1,0,0,0,1)},
        {-1,0,0,-2,0,0,2.001,4},
        {-1,1,0,-1,2,0,0,5},
        {-1,0,0,-3,0,0,2,4},
        {-1,-1,0,0,-1,0,0,5},
        {-1,-2,0,0,-2,0,0,5},
        {0,-1,0,1,0,0,0,5},
        {1,-1,0,1,0,0,0,5},
        {1,1,0,0,1,0,0,5},
        {0,1,0,0,2,0,0,5},
        {0,2,0,1,1,0,0,5},
        {0,-1,0,isr2,-1-isr2,0,(2-2*isr2)+0.001,4}
        };

        //std::vector<std::vector<float>> world={{-1,1,0, 1,1,0, 0,5}, {1,0,0, 2,0,0, 2,4}};

        //std::vector<std::vector<float>> world={{-1,1,0, 1,1,0, 0,5}};

        //std::vector<std::vector<float>> world={{2,0,0,4,0,0,0,mkdest(1,0,0,1,1)},{-1*isr2,1*isr2,0,-3*isr2,3*isr2,0,0,mkdest(1,0,1,1,1)},{0,-3,0,0,-2,0,3,4}};

        //std::vector<std::vector<float>> world={{4,0,0,2,0,0,2,mkdest(0,0,1,0,0)}, {-2,0,0,-4,0,0,2,mkdest(0,0,0,0,0)}};

        //std::vector<std::vector<float>> world={{3,0,0,1,0,0,2,mkdest(0,0,0,0,0)}};

        //std::vector<std::vector<float>> world={{0,0,0,6,0,0,4,mkdest(0,0,0,1,1)}};

        //std::vector<std::vector<float>> world={{3,0,0,1,0,0,2,mkdest(0,0,1,0,0)}, {-1,0,0,-3,0,0,2,mkdest(0,0,0,0,0)}
        //,{3,2,0,4,2,0,4,4}, {3,-2,0,4,-2,0,4,4}, {-1,2,0,-2,2,0,4,4}//, {-1,-2,0,-2,-2,0,4,4}
        //};
        //std::vector<std::vector<float>> world={{0,0,0,0,1,0,3,mkdest(0,0,1,0,0)},{3,0,0,2,0,0,3,mkdest(0,0,0,0,0)}};

        //std::vector<std::vector<float>> world={{2,1,0,0,1,0,3,mkdest(1,1,1,0,0)},{2,-1,0,0,-1,0,3,mkdest(1,1,0,0,1)},{-2,0,0,-3,0,0,2,4}};
        pawbuffer.push_back(world);
        worldCurvatures.push_back(0);
    }//*/
    /*{
        //std::vector<std::vector<float>> world={{2,1,0,0,1,0,3,mkdest(1,0,1,0,1)},{2,-1,0,0,-1,0,3,mkdest(1,0,0,0,0)}};
        std::vector<std::vector<float>> world={{0,1,0,-2,3,0,0,mkdest(1,0,2,0,0)},{-2,0,0,-1,0,0,4,mkdest(0,0,3,1,0)},{1,1,0,-1,-1,0,0,5},{1,1,0,2,0,0,0,5}};
        pawbuffer.push_back(world);
        worldCurvatures.push_back(0);
    }//*/

    //ms
    //00    side=1, mirror=0, invert location
    //10


    paw[0]=pawbuffer.size()+1;//serializing, redundant to later doing it every frame. This is so I can see debug info and be lazy
    for(int w=0;w<pawbuffer.size();w++){
        paw[w+1]=paw[w]+pawbuffer[w].size()*sizeOfPow;
        for(int i=0;i<pawbuffer[w].size();i++){
            if(pawbuffer[w][i].size() < 14) addExtraToPOW2(pawbuffer[w][i], worldCurvatures[w]);
            for(int p=0;p<sizeOfPow;p++){
                    paw[(int)paw[w]+sizeOfPow*i+p]=pawbuffer[w][i][p];
            }
        }
    }



    /*pl[0]=0;pl[1]=0;pl[2]=0;
    pl2[0]=pl[0]+pr;pl2[1]=pl[1];pl2[2]=pl[2]=0;
    plp[0]=1;
    camRef[0]=pl[0]+1;camRef[1]=pl[1];camRef[2]=pl[2];
    duppw = -1;
    pw=0;//*/
    //double testpt1[]={1,-1};
    //entity player;

    //pl[0]=0;pl[1]=0;pl[2]=1;pl2[0]=sin(pr);pl2[1]=0;pl2[2]=cos(pr);camRef[0]=1;camRef[1]=0;camRef[2]=0;
    //vec3(pl, 0,0,1); vec3(camRef, 1,0,sqrt(2));


//msiiiiiiiiiiwwwwwwwwwttt
//321098765432109876543210


    std::string shad="";
    {
        std::ifstream shaderCode("fs.glsl");
        char c;
        if(shaderCode.is_open()){
            while(shaderCode.get(c)) shad+=c;
        }
        shaderCode.close();
    }//*/

    /*std::ifstream shaderCode("fs.glsl");
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::string shad=buffer.str();*/


    std::string postprocessString="";
    {
        std::ifstream shaderCode("postProcess.glsl");
        char c;
        if(shaderCode.is_open()){
            while(shaderCode.get(c)) postprocessString+=c;
        }
        shaderCode.close();
    }

    char *fragmentShaderSource = &(shad[0]);
    char *postprocesscharpointer = &(postprocessString[0]);

    int width, height;

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    GLFWwindow* window = glfwCreateWindow(640, 480, "NightGun", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    glViewport(0, 0, 640, 480);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    float vertices[] = {
        1.0f, 1.0f, 0.0f,
        1.0f, -3.0f, 0.0f,
        -3.0f, 1.0f, 0.0f
    };
    unsigned int VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    unsigned int vertexShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    unsigned int fragmentShader;
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    GLint isCompiled = 0;
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &isCompiled);
    if(isCompiled == GL_FALSE){
        GLint maxLength = 0;
        glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<char> errorLog(maxLength);
        glGetShaderInfoLog(fragmentShader, maxLength, &maxLength, &errorLog[0]);

        //std::string* s=&errorLog[0];

        for(int n=0;n<maxLength;n++) std::cout<<errorLog[n];
        abort();
}


    unsigned int shaderProgram;
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);
    //
    unsigned int postprocesscode;
    postprocesscode = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(postprocesscode, 1, &postprocesscharpointer, NULL);
    glCompileShader(postprocesscode);
    //

    //
    unsigned int postprocessfull;
    postprocessfull = glCreateProgram();
    glAttachShader(postprocessfull, vertexShader);
    glAttachShader(postprocessfull, postprocesscode);
    glLinkProgram(postprocessfull);

    /*GLint program_linked;
    glGetProgramiv(postprocessfull, GL_LINK_STATUS, &program_linked);
    if (program_linked != GL_TRUE) printf("error\n");*/

    unsigned int framestage1;
    glGenTextures(1, &framestage1);
    //
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glDeleteShader(postprocesscode);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    unsigned int VAO;
    glGenVertexArrays(1, &VAO);
    // ..:: Initialization code (done once (unless your object frequently changes)) :: ..
    // 1. bind Vertex Array Object
    glBindVertexArray(VAO);
    // 2. copl[1] our vertices array in a buffer for OpenGL to use
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    // 3. then set our vertex attributes pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    int axesCount;
    int buttonCount;
    float deadzone = 0.06;
    float axi[6];
    GLuint ssbosd;
    glGenBuffers(1, &ssbosd);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbosd);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(shader_data), shader_data,  GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ssbosd);
    //glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // unbind
    unsigned int framebuffer1;
    glGenFramebuffers(1, &framebuffer1);
    glBindTexture(GL_TEXTURE_2D, framestage1);
    glUseProgram(shaderProgram);
    float zoom=1;
    long frameCount=0;
    float dx=0,dy=0;
    char oneshot=0;
    float facingAngle[5]={0,0,0,1,0};//3rd one is just working space



    double LOG[120][3+3+1+1 + 3+1+1 +1+1];//pl,cref, pw,iopip, duppl,duppw,diopip, dx,dy
    int loge = 0;
    while(!glfwWindowShouldClose(window)){
        double deltaTime=glfwGetTime();
        double NFT = deltaTime+1.0/60.0;//next frame time

        // input
        if(frameCount>0){
            double playerspeed = ps;
            bool running = false;
            int kx = 0, ky = 0;
            if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, true);
            if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) ky+=1;
            if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) kx-=1;
            if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) ky-=1;
            if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) kx+=1;//*/
            if(glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
                printP(pl);
                //printP(pl2);
                printP(camRef);
                //printP(duppl);
                //printf("%i\t%i\n",pw,duppw);
            }
            int logep = loge;
            if(glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) {
                loge--;if(loge < 0) loge = 120-1;
            }
            if(glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
                loge++;if(loge >= 120) loge = 0;
            }
            if(logep != loge){
                printf("%i\n",loge);
                for(int n = 0; n < 15; n++) printf("%lf\t",LOG[loge][n]);
                printf("\n");
                pl[0] = LOG[loge][0]; pl[1] = LOG[loge][1]; pl[2] = LOG[loge][2];
                camRef[0] = LOG[loge][3]; camRef[1] = LOG[loge][4]; camRef[2] = LOG[loge][5];
                pw = LOG[loge][6]; iopip = LOG[loge][7];
                duppl[0] = LOG[loge][8]; duppl[1] = LOG[loge][9]; duppl[2] = LOG[loge][10];
                duppw = LOG[loge][11];
                diopip = LOG[loge][12];
            }
            if(glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
                //-0.085489       0.996169        -0.018430       0.994265        0.086489        0.062896        1.000000        -1.000000       0.133523        0.496631        0.857630        -1.000000       -1.000000       0.000000        -0.032000
                vec3(pl, -0.085489,       0.996169,        -0.018430);
                vec3(camRef, 0.994265,        0.086489,        0.062896);
                pw = 1;
                moveEntity(0,-0.032000,pl,pl2,camRef,&pw,plp);
            }


            facingAngle[2]=0;
            const float *axes = glfwGetJoystickAxes(GLFW_JOYSTICK_1, &axesCount);
            //for(int n = 0; n < axesCount; n++) printf("%f\t",axes[n]);
            //printf("\n");
            if(axesCount > 0){
                if(abs(axes[0]) < deadzone) axi[0] = 0;
                else axi[0] = axes[0];
                if(abs(axes[1]) < deadzone) axi[1] = 0;
                else axi[1] = axes[1];
                if(abs(axes[rsxi]) < deadzone) axi[rsxi] = 0;
                else axi[rsxi] = axes[rsxi];
                if(abs(axes[rsyi]) < deadzone) axi[rsyi] = 0;
                else axi[rsyi] = axes[rsyi];
                dx = ps*axi[0];
                dy = -ps*axi[1];
                facingAngle[0] = axi[rsxi];
                facingAngle[1] = -axi[rsyi];
                if(abs(axi[rsxi]) > 0 || abs(axi[rsyi]) > 0) facingAngle[2] = 1;
                const unsigned char* buttons = glfwGetJoystickButtons(GLFW_JOYSTICK_1, &buttonCount);
                for(int n = 0; n < buttonCount; n++){
                    if(buttons[n] == GLFW_PRESS) printf("BUTTON: %i\n",n);
                }
                if(buttons[9] == GLFW_PRESS) running = true;
            }
            else {
                //for(int n = 0; n < 6; n++) axi[n] = 0;
                dx = dy = 0;
            }
            //A:0
            //B:1
            //X:2
            //Y:3
            //LB:4
            //RB:5
            //back:6
            //start:7
            //logo:8
            //LS:9
            //RS:10
            //D up:11
            //D right:12
            //D down:13
            //D left:14
            //

            if(kx != 0 || ky !=0){
                dx = kx*ps;
                dy = ky*ps;
                if(kx != 0 && ky != 0){
                    dx /= sr2;
                    dy /= sr2;
                }
            }

            if(facingAngle[2] > 0){
                //printf("%f\t%f\n",facingAngle[0],facingAngle[1]);
                facingAngle[2] = sqrt(facingAngle[0]*facingAngle[0] + facingAngle[1]*facingAngle[1]);
                if(facingAngle[2] > deadzone*2){
                    //printf("%f\n",facingAngle[2]);
                    facingAngle[3]=facingAngle[0]/facingAngle[2];
                    facingAngle[4]=facingAngle[1]/facingAngle[2]*plp[0];
                    //printf("%lf\t%lf\t%lf\n",facingAngle[2],facingAngle[3],facingAngle[4]);
                    if(worldCurvatures[pw]==0){//to do: do something similar for the duplicate
                        double crc[2]={camRef[0]-pl[0],camRef[1]-pl[1]};
                        pl2[0]=pr*facingAngle[3];pl2[1]=pr*facingAngle[4];
                        rotate(pl2[0],pl2[1],crc[0],crc[1]);
                        pl2[0]+=pl[0];pl2[1]+=pl[1];
                        pl2[2]=0;
                    } else if(worldCurvatures[pw]>0){
                        double rot[3][3],roti[3][3];
                        s2matto(pl,rot);
                        matxpt(rot,camRef,pl2);//reusing memory space that will be redefined anyway
                        rotXY(rot,pl2[0],-pl2[1]);
                        s2dadtopt(pr,facingAngle[3],facingAngle[4],pl2);
                        transpose(rot,roti);
                        matxpt(roti,pl2);
                        backOnSphere(pl2);
                    } else if(worldCurvatures[pw] < 0){
                        double rot[3][3],roti[3][3];
                        h2matto(pl,rot);
                        matxpt(rot,camRef,pl2);//reusing memory space that will be redefined anyway
                        rotXY(rot,pl2[0],-pl2[1]);
                        h2dadtopt(pr,facingAngle[3],facingAngle[4],pl2);
                        h2invert(rot,roti);
                        matxpt(roti,pl2);
                        backOnHyperboloid(pl2);
                    }
                    if(duppw >= 0 && dx == 0 && dy == 0){
                        updateDuplicateRot(pl,pl2,duppl,duppl2,pw,iopip);//TODO update this to work with H2
                        //banaid V
                        //moveEntity(ped,0,pl,pl2,camRef,&pw,plp);
                        //moveEntity(-ped,0,pl,pl2,camRef,&pw,plp);
                    }
                }
                facingAngle[2]=0;
            }
            if(dy!=0||dx!=0) {
                //movePlayer(dx,dy,0);
                if(running){
                    dx *= 2;
                    dy *= 2;
                }
                //pl,cref, pw,iopip, duppl,duppw,diopip
                loge++;if(loge >= 120) loge = 0;
                LOG[loge][0] = pl[0];LOG[loge][1] = pl[1];LOG[loge][2] = pl[2];
                LOG[loge][3] = camRef[0];LOG[loge][4] = camRef[1];LOG[loge][5] = camRef[2];
                LOG[loge][6] = pw; LOG[loge][7] = iopip;
                LOG[loge][8] = duppl[0]; LOG[loge][9] = duppl[1]; LOG[loge][10] = duppl[2];
                LOG[loge][11] = duppw;
                LOG[loge][12] = diopip;
                LOG[loge][13] = dx; LOG[loge][14] = dy;
                moveEntity(dx,dy,pl,pl2,camRef,&pw,plp);
                //printP(pl);
                //printP(camRef);
                //printP(duppl);
                //printf("%i\t%i\n",pw,duppw);
            }
            if(glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && zoom<1) zoom+=pow(2,floor(log2(zoom)-4));
            if(glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS && zoom>pr/pi) zoom-=pow(2,floor(log2(zoom)-4));
            glUniform3fv(glGetUniformLocation(shaderProgram,"pl"),1,&pl[0]);
            glUniform3fv(glGetUniformLocation(shaderProgram,"pl2"),1,&pl2[0]);
            glUniform3fv(glGetUniformLocation(shaderProgram,"camRef"),1,&camRef[0]);
            glUniform3fv(glGetUniformLocation(shaderProgram,"duppl"),1,&duppl[0]);
            glUniform3fv(glGetUniformLocation(shaderProgram,"duppl2"),1,&duppl2[0]);
            glUniform1f(glGetUniformLocation(shaderProgram, "mirrorPlayer"),plp[0]);
            glUniform1f(glGetUniformLocation(shaderProgram, "mirrorDup"),mirrorDup);
            glUniform1f(glGetUniformLocation(shaderProgram, "zoom"),zoom);
            glUniform1i(glGetUniformLocation(shaderProgram, "pw"),pw);
            glUniform1i(glGetUniformLocation(shaderProgram, "duppw"),duppw);
            glUniform1i(glGetUniformLocation(shaderProgram, "iopip"),iopip);
            glUniform1i(glGetUniformLocation(shaderProgram, "diopip"),diopip);
        }



        int widthp=width, heightp=height;
        glfwGetFramebufferSize(window, &width, &height);
        if(width!=widthp||height!=heightp) {
            glUniform1i(glGetUniformLocation(shaderProgram, "res"),std::min(width,height));
            printf("window dim:\n%i\tx\t%i\n",width,height);
            for(int n=0;n<paw[(int)paw[0]-1];n++) {
                if((int)paw[n]==paw[n]) printf("%i\t%i\n",n,(int)paw[n]);
                else printf("%i\t%f\n",n,paw[n]);
            }
        }
        // rendering commands here

        //
        //glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        //glClear(GL_COLOR_BUFFER_BIT);
        float T=1+sin(17.2787595947*glfwGetTime())/2;
        {
            for(int w=0;w<worldCurvatures.size();w++){
                worldCurvs[w]=worldCurvatures[w];
            }
            glUniform1fv(glGetUniformLocation(shaderProgram, "iaunits"),16,worldCurvs);
            shader_data[0]=pawbuffer.size()+1;//serializing
            for(int w=0;w<pawbuffer.size();w++){
                shader_data[w+1]=shader_data[w]+pawbuffer[w].size()*sizeOfPow;
                for(int i=0;i<pawbuffer[w].size();i++){
                    for(int p=0;p<sizeOfPow;p++){
                            shader_data[(int)shader_data[w]+sizeOfPow*i+p]=pawbuffer[w][i][p];
                    }
                }
            }
            shader_data[66559]=T;
            glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(shader_data), shader_data);
            //if(frameCount%60==0) std::cout<<glfwGetTime()-NFT+1/60.0<<"\n";
        }
        glUniform1f(glGetUniformLocation(shaderProgram, "T"),T);
        //int pawl = glGetUniformLocation(shaderProgram, "paw");
        //glUniform1fv(pawl,512,paw);
        //
        //1st pass
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framestage1, 0);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        //
        //2nd pass
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glUseProgram(postprocessfull);
        glUniform1i(glGetUniformLocation(postprocessfull, "height"),height);
        glUniform1i(glGetUniformLocation(postprocessfull, "width"),width);
        glUniform1f(glGetUniformLocation(postprocessfull, "zoom"),zoom);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        frameCount++;
        deltaTime=glfwGetTime()-deltaTime;

        // check and call events and swap the buffers
        glfwPollEvents();
        glfwSwapBuffers(window);
        glUseProgram(shaderProgram);

        while(glfwGetTime()<NFT);
        if(glfwGetTime()>=60) glfwSetTime(0);
    }

    glfwTerminate();
    //Close game controller
    //SDL_CloseJoystick( gGameController );
    //gGameController = NULL;

    //Quit SDL subsystems
    //SDL_Quit();
    return 0;
}

//to compile: g++ -O3 -w code.cpp -o code glad.o -lglfw -lSDL3 && ./code

//warnings are turned off by the -w, they were annoying

//helpful places:
//https://learnopengl.com/Getting-started/Hello-Triangle
//https://askubuntu.com/questions/1186517/which-package-to-install-to-get-header-file-glad-h
//https://registry.khronos.org/OpenGL-Refpages/gl4/html/glUniform.xhtml
//https://www.khronos.org/opengl/wiki/Shader_Storage_Buffer_Object
//https://stackoverflow.com/questions/71871742/how-can-i-get-joystick-data-using-glfw


//maybe https://gist.github.com/jasonwhite/c5b2048c15993d285130 ?


//NOT helpful places:
//StackOverflow
/*

"    int N=2,n=0;\n"
"    while(n<N&&zx*zx+zy*zy<4){\n"
"        t=zx*x-zy*y;zy=zx*y+zy*x;zx=t;\n"
"        n++;\n"
"    }\n"
"    i=float(n)/float(N);\n"
 */
