//to compile for windows on Linux: x86_64-w64-mingw32-g++ -O3 -w code.cpp -o nightgun.exe -static -L. -l glfw3dll
//to compile on Windows: g++ -O3 -w code.cpp -o nightgun.exe -static -L. -l glfw3dll
//^to get g++ working on Windows, follow this but only inclusively to "Check your MinGW installation": https://code.visualstudio.com/docs/cpp/config-mingw
//to run, open the folder this file is in, alt+d, type "cmd", hit enter, copy the compile command, and hit enter

//needed: sudo apt install libglfw3 libglfw3-dev

//for linux: g++ -O3 -w code.cpp -o nightgun -lglfw
#include <iostream>
#include "glad.h"
#include "glad.c"
#include "glfw3.h"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>

//#include <windows.h>


#define positionsSize 1048576

//SDL_Joystick* gGameController = NULL


/*std::vector<float> etfl(entity E){
    std::vector<float> F={E.pos1[0],E.pos1[1],E.pos1[2],E.pos2[0],E.pos2[1],E.pos2[2], E.props1, E.props2};
    return F;
}*/

std::vector<std::vector<std::vector<float>>> pawbuffer;
//std::vector<std::vector<std::vector<double>>> entitybuffer;//things that move around are stored with more precision
//std::vector<std::vector<std::vector<double>>> dopgangbuf;//doppleganger buffer, used for render and collision checking copies, not used yet
std::vector<float> worldCurvatures;
float shader_data[positionsSize];
float worldCurvs[128];//to do: make this all done in 1 batch later

const double tickRate = 60;
const float pi=3.14159265359, sr2=sqrt(2), isr2=sr2/2,ps = 0.043,pr=1.0/9.0, ped=1.0/1000, inf = 1/0.0;//player speed. player radius, portal ejection distance
const double wallThickness = 1.0/36;
const float prpwt=1.0/36+pr, prppwtc=cos(prpwt), prppwts=sin(prpwt);
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
double pl[3], pl2[3], camRef[3];
double playerduploc[3],playerduploc2[3];
float plp[1],paw[positionsSize];//player location, ..., player properties(mirror view flag (flip y)), duplicate player location, used for display purposes ,portals and walls
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

double e2side(double p1[], double p2[], double p3[]){
    return sign( (p3[0] - p1[0])*(p1[1] - p2[1]) + (p3[1] - p1[1])*(p2[0]) - p1[0]);
}


double disSquaredXY(double a[], double b[]){
    double c[2] = {a[0] - b[0], a[1] - b[1]};
    return c[0]*c[0] + c[1]*c[1];
}

double disSquared(double a[], double b[]){
    double c[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
}


void normalizeXY(double p[]){
    double r = sqrt(p[0]*p[0] + p[1]*p[1]);
    p[0] /= r;
    p[1] /= r;
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

void copypt(float p1[3], float p2[3]){
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

void matxpt(double m[][3], float p[3], double p2[3]){
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

void rotateXY(double m[][3], double c, double s){//to do, swap order of inputs
    rotate(m[0][0],m[0][1],c,s);
    rotate(m[1][0],m[1][1],c,s);
    rotate(m[2][0],m[2][1],c,s);
}

void rotateXZ(double m[][3], double c, double s){
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
    rotateXY(m,c2,-s2);
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

void lorentz(double& x, double& y, double c, double s){
    double t = x*c + y*s;
    y = x*s + y*c;
    x=t;
}

void lorentz(float& x, float& y, double c, double s){
    double t = x*c + y*s;
    y = x*s + y*c;
    x=t;
}

void lorentzXZ(double m[][3], double c, double s){
    lorentz(m[0][0],m[0][2],c,s);
    lorentz(m[1][0],m[1][2],c,s);
    lorentz(m[2][0],m[2][2],c,s);
}

void lorentzYZ(double m[][3], double c, double s){
    lorentz(m[0][1],m[0][2],c,s);
    lorentz(m[1][1],m[1][2],c,s);
    lorentz(m[2][1],m[2][2],c,s);
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

void idealrot2(double m[3][3], double d){//TODO make invert the signs of all d's
    double A = d*d/2, B = 1 - A, C = 1 + A, a, b, c;
    for(int n = 0; n < 3; n++){
        a = m[n][0]*B - m[n][1]*d - m[n][2]*A;
        b = m[n][1] + (m[n][0] + m[n][2])*d;
        c = m[n][2]*C + m[n][0]*A + m[n][1]*d;
        vec3(m[n],a,b,c);
    }
    //...
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

void backOnHyperboloid(double p[]){
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

double disSquaredAprox(double a[], double b[], double iaunit){
    double d2;
    if(iaunit == 0){
        double c[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
        d2 = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
    } else if(iaunit > 0){
        d2 = 2 - 2*dot(a,b);
    } else if(iaunit < 0) {
        d2 = 2*lidot(a,b) - 2;
    }
    return d2;
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

std::vector<float> s2threepointstowall(double p5[], double p3[], double p4[]){
    //getting normal
    double p1[3], p3c[3], p4c[3];
    vec3(p3c, p3[0] - p5[0],p3[1] - p5[1],p3[2] - p5[2]);
    vec3(p4c, p4[0] - p5[0],p4[1] - p5[1],p4[2] - p5[2]);
    cross(p3c,p4c,p1);
    backOnSphere(p1);
    //getting mitpoint
    double m1[3][3], m2[3][3];
    s2matto(p1,m1);
    matxpt(m1,p3);
    matxpt(m1,p4);
    matxpt(m1,p5);
    double r = sqrt(1 - p3[2]*p3[2]);
    double p2[] = {p3[1] - p4[1],p4[0] - p3[0],p4[2]};
    double s = r/sqrt(p2[0]*p2[0] + p2[1]*p2[1]);
    p2[0] *= s; p2[1] *= s;
    //correcting if wrong midpoint;
    p3[0] -= p2[0]; p3[1] -= p2[1];
    p5[0] -= p2[0]; p5[1] -= p2[1];
    if(p5[0]*p5[0] + p5[1]*p5[1] > p3[0]*p3[0] + p3[1]*p3[1]){
        p2[0] *= -1; p2[1] *= -1;
    }
    double limit = (dot(p2,p4) - 1.00001)/r/r;
    transpose(m1,m2);
    matxpt(m2,p2);
    backOnSphere(p2);
    std::vector<float> wall = {p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], limit,4};
    return wall;
    //...
}

std::vector<float> h2threepointstowall(double p2[], double p3[], double p4[]){
    double c,s,d;
    double p3c[3], p4c[3];
    copypt(p3,p3c);
    copypt(p4,p4c);
    p3c[0] -= p2[0];p3c[1] -= p2[1];p3c[2] -= p2[2];
    p4c[0] -= p2[0];p4c[1] -= p2[1];p4c[2] -= p2[2];
    double p1[3];
    cross(p3c,p4c,p1);
    p1[2] *= -1;
    s = lidot(p1,p1);
    int m = sign(s);
    if(s != 0){
        s = sqrt(m/s);
    }
    p1[0] *= s;p1[1] *= s;p1[2] *= s;
    double k = lidot(p1,p2);
    double p1c[3];
    copypt(p1,p1c);
    copypt(p4,p4c);
    double m1[3][3],m2[3][3];
    h2matto(p3,m1);
    matxpt(m1,p4c);
    d = sqrt(p4c[0]*p4c[0]+p4c[1]*p4c[1]);
    c = p4c[0]/d;s = p4c[1]/d;
    rotateXY(m1,-s,-c);
    rotate(p4c[0],p4c[1],-s,-c);
    c = lidot(p3,p4);
    c = sqrt((c+1)/2);
    s = sqrt(c*c-1);
    lorentzYZ(m1,c,s);
    matxpt(m1,p2);
    matxpt(m1,p1c);
    double B = k*p1c[0];
    p2[0] = (B+sign(p2[0])*sqrt(B*B+m*(k*k - p1c[2]*p1c[2])))/m;
    p2[1] = 0;
    p2[2] = sqrt(p2[0]*p2[0]+1);
    h2invert(m1,m2);
    matxpt(m2,p2);
    c = lidot(p2,p3);
    std::vector<float> wall = {p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], 1.0000001*c,5-m};
    return wall;
}

void addExtraToPOW(std::vector<float>& G, float iaunit){
    double p1[3] = {G[0],G[1],G[2]}, p2[3] = {G[3],G[4],G[5]};
    double L = G[6];
    int dat = (int)G[7];
    char type = dat&7;
    bool iswall = (type > 3);
    if(iaunit > 0){
        if(iswall || true){
            //
            double A,B,D,cr,sr2,t;
            cr = dot(p1,p2);
            sr2 = 1 - cr*cr;
            double L2 = 1 + sr2*L;//TODO remove the conversion and just store as the radius dependent value
            B = (L2 - cr*cr)/sr2;
            A = cr*(1 - B);
            D = A*A + 2*A*B*cr+B*B;
            t = sqrt( (1 - D)/sr2);
            double p3[3];
            cross(p1,p2,p3);
            double endpt1[3], endpt2[3];
            for(int n = 0; n < 3; n++) {
                endpt1[n] = A*p1[n] + B*p2[n] + t*p3[n];
                endpt2[n] = A*p1[n] + B*p2[n] - t*p3[n];
            }
            //*/
            /*/
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
            matxpt(roti,endpt2);//*/
            G.push_back(endpt1[0]);
            G.push_back(endpt1[1]);
            G.push_back(endpt1[2]);
            G.push_back(endpt2[0]);
            G.push_back(endpt2[1]);
            G.push_back(endpt2[2]);
            //printP(endpt1);
            //printP(endpt2);
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
        if(true){
            double k = lidot(p1,p2);
            double rot[3][3], roti[3][3];
            h2matto(p2,rot);
            matxpt(rot,p1);
            double r = sqrt(p1[0]*p1[0] + p1[1]*p1[1]);
            rotateXY(rot, p1[0]/r,-p1[1]/r);
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

void getTouchingChunksNary(){
    //
}


class entity{
public:
    int currentWorld;
    double cord1[3];
    double cord2[3];
    double compass[3];//expected to have length 1
    bool bounce;
    bool isDupe;
    double radius;
    float mirror;//could be a bool, but would cause more coversions
    //...
    int intersectingPortalIndex;
    entity* duplicate;
    //...
};

void printE(entity E){
    printf("cw:\t%i\n",E.currentWorld);
    printP(E.cord1);
    printP(E.cord2);
    printP(E.compass);
}

void copyEntityState(entity& E, entity& D){
    D.currentWorld = E.currentWorld;
    copypt(E.cord1,D.cord1);
    copypt(E.cord2,D.cord2);
    copypt(E.compass,D.compass);
    D.bounce = E.bounce;
    D.isDupe = E.isDupe;
    D.radius = E.radius;
    D.mirror = E.mirror;
    //...
    D.intersectingPortalIndex = E.intersectingPortalIndex;

    //D.duplicate = E.duplicate;
}

void clearDuplicates(entity& E){
    E.intersectingPortalIndex = -1;
    E.duplicate[0].intersectingPortalIndex = -1;
}

void addDuplicate(entity& E, int i, double relLoc, double disToCp, double rp[3][3]){
    //printf("%lf\t%lf\n",relLoc,disToCp);
    int dat = paw[i + 7];
    E.intersectingPortalIndex = i;
    E.duplicate[0].currentWorld = ((dat)>>3)&511;
    E.duplicate[0].radius = E.radius;
    E.duplicate[0].isDupe = true;
    int di = int(paw[E.duplicate[0].currentWorld]) + sizeOfPow*((dat>>12)&1023);
    E.duplicate[0].intersectingPortalIndex = di;
    double mirror, side;
    if(((dat>>22)&1)==1) side = -1;
    else side = 1;
    E.duplicate[0].mirror = E.mirror;
    if(((dat>>23)&1)==1) {mirror = -1; E.duplicate[0].mirror *= -1;}
    else mirror = 1;
    relLoc *= mirror*side;
    disToCp *= side;
    double r;
    for(int n = 0; n < 3; n++){
        r = sqrt(rp[n][0]*rp[n][0] + rp[n][1]*rp[n][1]);
        rp[n][0] /= r; rp[n][1] /= r;
    }
    double* p;
    rotate(rp[0][0],rp[0][1],-rp[2][0],rp[2][1]);
    rotate(rp[1][0],rp[1][1],-rp[2][0],rp[2][1]);
    double signdtcp = sign(disToCp);
    rp[0][0] *= signdtcp;
    rp[0][1] *= mirror*signdtcp;
    rp[1][0] *= signdtcp;
    rp[1][1] *= mirror*signdtcp;
    double iaunit = worldCurvatures[E.duplicate[0].currentWorld];
    dat = paw[di + 7];
    char type = char(dat&3);
    iaunit = worldCurvatures[E.duplicate[0].currentWorld];
    double sourcept1[] = {paw[di],paw[di+1],paw[di+2]};
    double sourcept2[] = {paw[di+3],paw[di+4],paw[di+5]};
    if(iaunit < 0){
        double k = lidot(sourcept1,sourcept2);
        if(type > 0){
            double tempvec[3];
            copypt(sourcept1,tempvec);
            copypt(sourcept2,sourcept1);
            copypt(tempvec,sourcept2);
        }
        double aunit = -1/iaunit;
        double transform[3][3], itransform[3][3];
        h2matto(sourcept1,transform);
        matxpt(transform,sourcept2);
        double r = sqrt(k*k - 1 + type);
        if(type > 0){
            sourcept2[0] *= -1;
            sourcept2[1] *= -1;
        }
        rotateXY(transform,sourcept2[0]/r,-sourcept2[1]/r);
        relLoc /= -aunit*r;
        if(type == 0){
            rotateXY(transform,cos(relLoc),sin(relLoc));
            lorentzXZ(transform,k,-r);
        } else if(type == 1){
            idealrot2(transform,relLoc*r);
        } else if(type == 2){
            double s = sinh(relLoc), c = sqrt(1 + s*s);
            lorentzXZ(transform,r,k);
            lorentzYZ(transform,c,s);
            lorentzXZ(transform,r,-k);
        }
        double s = sinh(disToCp/aunit), c = sqrt(1 + s*s);
        lorentzXZ(transform,c,-s);
        vec3(E.duplicate[0].cord1, 0,0,1);
        h2dadtopt(E.radius/aunit, rp[0][0],rp[0][1], E.duplicate[0].cord2);
        vec3(E.duplicate[0].compass, rp[1][0],rp[1][1],sr2);
        h2invert(transform,itransform);
        matxpt(itransform,E.duplicate[0].cord1);
        matxpt(itransform,E.duplicate[0].cord2);
        matxpt(itransform,E.duplicate[0].compass);
    }
    if(iaunit == 0){
        double off[2] = {-sourcept1[0],-sourcept1[1]};
        sourcept2[0] += off[0]; sourcept2[1] += off[1];
        r = sqrt(sourcept2[0]*sourcept2[0] + sourcept2[1]*sourcept2[1]);
        double rot[2] = {sourcept2[0]/r,-sourcept2[1]/r};
        rotate(off[0],off[1],rot[0],rot[1]);
        if(type == 0){
            double a = relLoc/r;
            double c = cos(a), s = - sin(a);
            rotate(rot[0],rot[1],c,s);
            rotate(off[0],off[1],c,s);
            off[0] -= r;
        } else if(type == 1){
            off[0] -= r/2 + relLoc;
            double t;
            t = rot[0]; rot[0] = -rot[1]; rot[1] = t;
            t = off[0]; off[0] = -off[1]; off[1] = t;
            //rotate the world 90 degrees^
        }
        off[0] -= disToCp;
        vec3(E.duplicate[0].cord1, 0,0,0);
        vec3(E.duplicate[0].cord2, E.radius*rp[0][0],E.radius*rp[0][1],0);
        vec3(E.duplicate[0].compass, rp[1][0],rp[1][1],0);
        p = E.duplicate[0].cord1;
        p[0] -= off[0]; p[1] -= off[1];
        rotate(p[0],p[1],rot[0],-rot[1]);
        p = E.duplicate[0].cord2;
        p[0] -= off[0]; p[1] -= off[1];
        rotate(p[0],p[1],rot[0],-rot[1]);
        p = E.duplicate[0].compass;
        p[0] -= off[0]; p[1] -= off[1];
        rotate(p[0],p[1],rot[0],-rot[1]);
    }
    if(iaunit > 0){
        double rotation[3][3];
        double cr = dot(sourcept1,sourcept2), sr = sqrt(1 - cr*cr);
        s2matto(sourcept1,rotation);
        matxpt(rotation,sourcept2);
        sourcept2[0] /= sr; sourcept2[1] /= sr;
        rotateXY(rotation, sourcept2[0],-sourcept2[1]);
        double a = -relLoc*iaunit/sr;
        rotateXY(rotation, cos(a),sin(a));
        rotateXZ(rotation, cr,sr);
        disToCp *= iaunit;
        rotateXZ(rotation,cos(disToCp),sin(disToCp));
        vec3(E.duplicate[0].cord1, 0,0,1);
        s2dadtopt(E.radius*iaunit, rp[0][0],rp[0][1], E.duplicate[0].cord2);
        vec3(E.duplicate[0].compass, rp[1][0],rp[1][1], 0);
        double irotation[3][3];
        transpose(rotation,irotation);
        matxpt(irotation,E.duplicate[0].cord1);
        matxpt(irotation,E.duplicate[0].cord2);
        matxpt(irotation,E.duplicate[0].compass);
    }
}

void updateDuplicateRot(entity& E, double joyX, double joyY){
    //...
    joyX *= E.radius;
    joyY *= E.radius;
    double iaunit = worldCurvatures[E.duplicate[0].currentWorld];
        //...
        if(iaunit < 0){
            joyX = -sinh(joyX*iaunit);
            joyY = -sinh(joyY*iaunit);
            double transform[3][3], itransform[3][3];
            h2matto(E.duplicate[0].cord1,transform);
            matxpt(transform,E.duplicate[0].compass,E.duplicate[0].cord2);
            rotate(E.duplicate[0].cord2[0],E.duplicate[0].cord2[1],joyX,joyY*E.duplicate[0].mirror*E.mirror);
            backOnHyperboloid(E.duplicate[0].cord2);
            h2invert(transform,itransform);
            matxpt(itransform,E.duplicate[0].cord2);
        }
        if(iaunit == 0){
            E.duplicate[0].cord2[0] = E.duplicate[0].compass[0] - E.duplicate[0].cord1[0];
            E.duplicate[0].cord2[1] = E.duplicate[0].compass[1] - E.duplicate[0].cord1[1];
            rotate(E.duplicate[0].cord2[0],E.duplicate[0].cord2[1],joyX,joyY*E.duplicate[0].mirror*E.mirror);
            E.duplicate[0].cord2[0] = E.duplicate[0].cord2[0] + E.duplicate[0].cord1[0];
            E.duplicate[0].cord2[1] = E.duplicate[0].cord2[1] + E.duplicate[0].cord1[1];
        }
        //...
    //...
}

void moveEntity(double dx, double dy, entity& E, bool& abortMove){
    const int maxits = 5;
    double d = sqrt(dx*dx + dy*dy);
    dx /= d; dy /= d*E.mirror;
    double mirror = 1;
    double side = 1;
    double entranceAng[2] = {dx,-dy};
    double entranceLoc;
    int type = 0;
    int sourceIndex = -1;
    double sourcept1[3], sourcept2[3],cord2Rel[3];
    copypt(E.cord1,sourcept1);copypt(E.compass,sourcept2);
    int iterations = 0;
    double transform[3][3];
    double rot[2], off[2];
    while(d > 0 && iterations < maxits && !abortMove){
        double iaunit = worldCurvatures[E.currentWorld];
        iterations++;
        int cii = -1;
        int di, dat;
        double cip1[3];
        if(iaunit < 0){//hyperbolic
            double k = lidot(sourcept1,sourcept2);
            double r = sqrt(k*k - 1 + type);
            double aunit = -1/iaunit;
            double ds = sinh(d/aunit);
            double dc = sqrt(ds*ds + 1);
            if(type < 4){
                if(type > 0){
                    double tempvec[3];
                    copypt(sourcept1,tempvec);
                    copypt(sourcept2,sourcept1);
                    copypt(tempvec,sourcept2);
                }
                h2matto(sourcept1,transform);
                matxpt(transform,sourcept2);
                if(type > 0){
                    sourcept2[0] *= -1;
                    sourcept2[1] *= -1;
                }
                rotateXY(transform,sourcept2[0]/r,-sourcept2[1]/r);
                entranceLoc /= -aunit*r;
                if(type == 0 && sourceIndex > 0){
                    rotateXY(transform,cos(entranceLoc),sin(entranceLoc));
                    lorentzXZ(transform,k,-r);
                } else if(type == 1){
                    idealrot2(transform,entranceLoc*r);
                    //...
                } else if(type == 2){
                    double s = sinh(entranceLoc), c = sqrt(1 + s*s);
                    lorentzXZ(transform,r,k);
                    lorentzYZ(transform,c,s);
                    lorentzXZ(transform,r,-k);
                }
                rotateXY(transform,entranceAng[0],entranceAng[1]);
                if(sourceIndex == -1){
                    matxpt(transform,E.cord2,cord2Rel);
                    normalizeXY(cord2Rel);
                }
            }
            //
            double cip2[3], cip[3];
            cip[0] = ds;
            int i = int(paw[E.currentWorld]);
            while(i < (int)paw[E.currentWorld + 1]){
                double p1[3] = {paw[i],paw[i + 1],paw[i + 2]},
                p2[3] = {paw[i + 3],paw[i + 4],paw[i + 5]};
                double limit = paw[i + 6];
                type = int(paw[i + 7])&3;
                k = lidot(p1,p2);
                matxpt(transform,p1);
                matxpt(transform,p2);
                double o = sqrt(k*k - p1[2]*p1[2] + p1[0]*p1[0]);
                if(p1[2] < 0) o = -o;
                o += k;
                double tL = (p1[2] - p1[0])/o;
                double tR = o/(p1[2] + p1[0]);
                double cipR[3], cipL[3];
                cipR[1] = cipL[1] = 0;
                cipR[0] = (1/tR - tR)/2; cipR[2] = cipR[0] + tR;
                cipL[0] = (1/tL - tL)/2; cipL[2] = cipL[0] + tL;
                bool risvalid = tR > 0 && limit > lidot(p2,cipR);
                bool lisvalid = tL > 0 && limit > lidot(p2,cipL);
                if(risvalid && 0 < cipR[0] && cipR[0] < cip[0] && (i != sourceIndex || (lisvalid && cipL[2] < cipR[2]))){
                    copypt(cipR,cip);
                    cii = i;
                }
                if(lisvalid && 0 < cipL[0] && cipL[0] < cip[0] && (i != sourceIndex || (risvalid && cipR[2] < cipL[2]))){
                    copypt(cipL,cip);
                    cii = i;
                }
                if(i == cii){
                    copypt(p1,cip1);
                    copypt(p2,cip2);
                }
                i += sizeOfPow;
            }
            //
            if(cii > 0){
                dat = int(paw[cii + 7]);
                type = dat&7;
                if(type > 3){
                    lorentzXZ(transform, cip[2],-cip[0]);
                    lorentz(cip1[0],cip1[2],cip[2],-cip[0]);
                    //
                    double r2 = cip1[0]*cip1[0] + cip1[1]*cip1[1];
                    double c = 2*cip1[1]*cip1[1]/r2 - 1, s = -2*cip1[0]*cip1[1]/r2;
                    //
                    rotate(cord2Rel[0],cord2Rel[1],c,-s);
                    rotateXY(transform,c,-s);
                    d -= aunit*acosh(cip[2]) - ped;
                    sourceIndex = cii;
                    //d = 0;
                } else {
                    double topbot = h2side(cip2,cip1,cip);
                    entranceLoc = lidot(cip,cip2);
                    if(entranceLoc < 1) entranceLoc = 1;
                    double k = lidot(cip1,cip2);
                    E.currentWorld = ((dat)>>3)&511;
                    di = int(paw[E.currentWorld]) + sizeOfPow*((dat>>12)&1023);
                    sourceIndex = di;
                    vec3(sourcept1, paw[di],paw[di + 1],paw[di + 2]);
                    vec3(sourcept2, paw[di + 3],paw[di + 4],paw[di + 5]);
                    if(((dat>>22)&1)==1) side=-1;
                    else side=1;
                    if(((dat>>23)&1)==1) {mirror = -1; E.mirror *=-1;}
                    else mirror=1;
                    r = k*k - 1 + type;
                    lorentz(cip1[0],cip1[2],cip[2],-cip[0]);
                    if(type == 0){
                        entranceLoc = acos((k*k - entranceLoc)/r);
                    } else if(type == 1){
                        entranceLoc = sqrt(entranceLoc*2 - 2)/k;
                    } else {
                        entranceLoc = acosh((k*k + entranceLoc)/r);
                    }
                    r = safe_sqrt(r)*side;
                    entranceLoc *= r*aunit*topbot*mirror;
                    entranceAng[0] = -cip1[0]/r; entranceAng[1] = -cip1[1]/r*mirror;
                    type = char(int(paw[di+7])&3);
                    d -= aunit*acosh(cip[2]) - ped;
                }
                //...
                /* Common code:
                E.currentWorld = ((dat)>>3)&511;
                di = int(paw[E.currentWorld]) + sizeOfPow*((dat>>12)&1023);
                vec3(sourcept1, paw[di],paw[di + 1],paw[di + 2]);
                vec3(sourcept2, paw[di + 3],paw[di + 4],paw[di + 5]);
                sourceIndex = di;
                if(((dat>>22)&1)==1) side=-1;
                else side=1;
                if(((dat>>23)&1)==1) {mirror = -1; E.mirror *=-1;}
                else mirror=1;
                r *= side;
                entranceAng[0] = -cip1[0]/r; entranceAng[1] = -cip1[1]/r*mirror;
                entranceLoc *= mirror*side;
                type = char(int(paw[di+7])&3);
                 */
            } else {
                lorentzXZ(transform,dc,-ds);
                double itransform[3][3];
                h2invert(transform,itransform);
                vec3(E.cord1, 0,0,1);
                vec3(E.compass,dx,-dy*mirror,sr2);
                h2dadtopt(E.radius/aunit, cord2Rel[0],cord2Rel[1]*mirror, E.cord2);
                matxpt(itransform,E.cord1);
                matxpt(itransform,E.cord2);
                matxpt(itransform,E.compass);
                backOnHyperboloid(E.cord1);
                backOnHyperboloid(E.cord2);
                backOnHyperboloid(E.compass);
                d = 0;
                //
                clearDuplicates(E);
                i=(int)paw[E.currentWorld];
                bool inWall = false;
                double cp[3];
                double coshEr = cosh(E.radius*iaunit);
                double coshEcolr = cosh((E.radius + wallThickness)*iaunit);
                double coshped = cosh(ped*iaunit);
                while(i < paw[E.currentWorld + 1] && !inWall){
                    type = int(paw[i + 7])&7;
                    bool isWall;
                    if(type > 3) isWall = true;
                    else isWall = false;
                    if(!(isWall && E.bounce)){
                        type &= 3;
                        double p1[] = {paw[i],paw[i + 1],paw[i + 2]}, p2[] = {paw[i + 3],paw[i + 4],paw[i + 5]};
                        double p0[] = {E.cord1[0],E.cord1[1],E.cord1[2]};
                        k = lidot(p1,p2);
                        double i01 = lidot(p0,p1), i11 = 1 - type;
                        double C = i01*i01 - i11;
                        double A = sqrt( (k*k - i11)/C );
                        double B;
                        if(type == 1) B = (i01*i01 - k*k)/(2*C*k);
                        else B = (k - A*i01)/i11;
                        for(int n = 0; n < 3; n++) cp[n] = A*p0[n] + B*p1[n];
                        double updown = h2side(p1,p2,p0);
                        double relLoc = lidot(cp,p2);
                        double limit = paw[i + 6];
                        if(relLoc > limit){
                            relLoc = limit;
                            int off;
                            if(updown > 0) off = 8;
                            else off = 11;
                            vec3(cp, paw[i + off],paw[i + off + 1],paw[i + off + 2]);
                        }
                        double coshDisToCp = lidot(cp,p0);
                        if( (isWall && coshDisToCp < coshEcolr) || (!isWall && coshDisToCp < coshped)){
                            double transform[3][3];
                            h2matto(p0,transform);
                            matxpt(transform,cp);
                            double compasscopy[3];
                            matxpt(transform,E.compass,compasscopy);
                            rotate(cp[0],cp[1],-compasscopy[0],compasscopy[1]);
                            normalizeXY(cp);
                            d = ped - acosh(coshDisToCp)*aunit;
                            if(isWall) d += E.radius + wallThickness;
                            else d += ped;
                            dx = cp[0]; dy = cp[1];
                            entranceAng[0] = dx; entranceAng[1] = -dy;
                            type = 0;
                            sourceIndex = -1;
                            copypt(E.cord1,sourcept1);
                            copypt(E.compass,sourcept2);
                            inWall = true;
                            mirror = 1;
                        }
                        r = k*k - 1 + type;
                        if(!inWall && coshDisToCp < coshEr && !isWall){
                            double relLoc = lidot(cp,p2);
                            if(type == 0){
                                relLoc = (k*k - relLoc)/r;
                                if(relLoc > 1) relLoc = 1;
                                if(relLoc < -1) relLoc = -1;
                                relLoc = acos(relLoc);
                            } else if(type == 1){
                                relLoc = safe_sqrt(2*relLoc - 2)/k;
                            } else {
                                relLoc = (k*k + relLoc)/r;
                                if(relLoc < 1) relLoc = 1;
                                relLoc = acosh(relLoc);
                            }
                            r = sqrt(r);
                            relLoc *= -aunit*r*updown;
                            double transform[3][3];
                            h2matto(p0,transform);
                            double relPoses[3][3];
                            copypt(E.cord2,relPoses[0]);
                            copypt(E.compass,relPoses[1]);
                            copypt(cp,relPoses[2]);
                            for(int n = 0; n < 3; n++) matxpt(transform,relPoses[n]);
                            double disToCp = acosh(coshDisToCp)*sign(lidot(p0,p1)-k)*aunit;
                            addDuplicate(E,i,relLoc,disToCp,relPoses);
                        }
                    }
                    i += sizeOfPow;
                }
                //*/
            }
            //...

        } else if(iaunit == 0){//euclidian
            double r;
            //
            if(type < 4){
                off[0] = -sourcept1[0]; off[1] = -sourcept1[1];
                sourcept2[0] += off[0]; sourcept2[1] += off[1];
                r = sqrt(sourcept2[0]*sourcept2[0] + sourcept2[1]*sourcept2[1]);
                rot[0] = sourcept2[0]/r; rot[1] = -sourcept2[1]/r;
                rotate(off[0],off[1],rot[0],rot[1]);
                if(type == 0 && sourceIndex > 0){
                    double a = entranceLoc/r;
                    double c = cos(a), s = - sin(a);
                    rotate(rot[0],rot[1],c,s);
                    rotate(off[0],off[1],c,s);
                    off[0] -= r;
                } else if(type == 1){
                    off[0] -= r/2 + entranceLoc;
                    double t;
                    t = rot[0]; rot[0] = -rot[1]; rot[1] = t;
                    t = off[0]; off[0] = -off[1]; off[1] = t;
                    //rotate the world 90 degrees^
                }
                rotate(rot[0],rot[1],entranceAng[0],entranceAng[1]);
                rotate(off[0],off[1],entranceAng[0],entranceAng[1]);
                if(sourceIndex == -1){
                    cord2Rel[0] = E.cord2[0]; cord2Rel[1] = E.cord2[1];
                    rotate(cord2Rel[0],cord2Rel[1],rot[0],rot[1]);
                    cord2Rel[0] += off[0]; cord2Rel[1] += off[1];
                    normalizeXY(cord2Rel);
                }
            }
            int i=(int)paw[E.currentWorld],cii=-1;
            double cip2[2];
            double cid=d;
            char cit;
            while(i < paw[E.currentWorld + 1]){
                type = ((char)paw[i+7])&7;
                char typetemp = type;
                if(type < 4 || E.bounce){
                    type &= 3;
                    double p1[] = {paw[i],paw[i + 1]}, p2[] = {paw[i + 3],paw[i + 4]};
                    r = (p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]);//actually r^2
                    double al = paw[i + 6];//TODO store things with al already multiplied by r^2
                    rotate(p1[0],p1[1],rot[0],rot[1]);
                    p1[0] += off[0];p1[1] += off[1];
                    rotate(p2[0],p2[1],rot[0],rot[1]);
                    p2[0] += off[0];p2[1]+=off[1];
                    if(type == 0){
                        if(p1[1]*p1[1]<=r&&(al==4*r||p2[1]*p2[1]<=al*r)){
                            double o=sqrt(r-p1[1]*p1[1]), xt=p1[0]-o, xt1=abs(xt);
                            if(((al==4) || ((xt-p2[0])*(xt-p2[0])+p2[1]*p2[1]<=al*r)) && xt>0 && xt<cid && i!=sourceIndex ){
                                cip1[0]=p1[0];cip1[1]=p1[1];
                                cip2[0]=p2[0];cip2[1]=p2[1];
                                cii = i; cid = xt;
                            }
                            xt=p1[0]+o;
                            if(xt>0){
                                if(((al==4) || ((xt-p2[0])*(xt-p2[0])+p2[1]*p2[1]<=al*r)) && xt<cid && (i!=sourceIndex || xt>xt1) ){//problem place, think it is fixed now.
                                    cip1[0]=p1[0];cip1[1]=p1[1];
                                    cip2[0]=p2[0];cip2[1]=p2[1];
                                    cii = i; cid = xt;
                                }
                            }
                        }
                    } else {
                        if((p1[1]*p2[1]<0)&&((p1[0]>=0)||(p2[0]>=0))&&(i!=sourceIndex)){//crosses x axis and one of the 2 points is in front and not the source portal
                            double xt=(p1[0]-p2[0])*p2[1]/(p2[1]-p1[1])+p2[0];
                            if(xt>0 && xt<cid){
                                cip1[0]=p1[0];cip1[1]=p1[1];
                                cip2[0]=p2[0];cip2[1]=p2[1];
                                cii = i; cid = xt;
                            }
                        }
                    }
                    if(cii == i) cit = typetemp;
                }
                i+=sizeOfPow;
            }
            if(cii > 0){//hit a wall or portal
                //
                if(cit < 4){
                    cip1[0] -= cid; cip2[0] -= cid;
                    r = -sqrt(cip1[0]*cip1[0] + cip1[1]*cip1[1]);
                    entranceAng[0] = cip1[0]/r; entranceAng[1] = cip1[1]/r;
                    cip2[0] -= cip1[0]; cip2[1] -= cip1[1];
                    cip1[0] *= -1; cip1[1] *= -1;
                    double R2 = cip2[0]*cip2[0] + cip2[1]*cip2[1];
                    double R = sqrt(R2);
                    if(cit == 0){
                        entranceLoc = R*sign(cip1[1]*cip2[0] - cip1[0]*cip2[1])*acos(1 - disSquaredXY(cip1,cip2)/2/R2);
                    } else {
                        double t = entranceAng[0];
                        entranceAng[0] = entranceAng[1];
                        entranceAng[1] = -t;
                        entranceLoc = -r- R/2;
                    }
                    dat = (int)paw[cii + 7];
                    E.currentWorld = ((dat)>>3)&511;
                    di = int(paw[E.currentWorld]) + sizeOfPow*((dat>>12)&1023);
                    sourceIndex = di;
                    sourcept1[0]=paw[di];sourcept1[1]=paw[di+1];sourcept1[2]=paw[di+2];
                    sourcept2[0]=paw[di+3];sourcept2[1]=paw[di+4];sourcept2[2]=paw[di+5];
                    type = char(int(paw[di+7])&7);
                    if(((dat>>22)&1)==1) side=-1;
                    else side=1;
                    if(((dat>>23)&1)==1) {mirror = -1; E.mirror *=-1;}
                    else mirror=1;
                    entranceAng[0] *= side; entranceAng[1] *= side*mirror;
                    entranceLoc *= mirror*side;
                    d -= cid - ped;
                } else {
                    type = cit;
                    off[0] -= cid;
                    cip1[0] -= cid;
                    if(cit == 5){
                        cip1[2] = cip1[0];//unused in E2, so using it as scratch memory
                        cip1[0] = -cip1[1];
                        cip1[1] = cip1[2];
                    }
                    //
                    double r2 = cip1[0]*cip1[0] + cip1[1]*cip1[1];
                    double c = 2*cip1[1]*cip1[1]/r2 - 1, s = -2*cip1[0]*cip1[1]/r2;
                    //
                    rotate(cord2Rel[0],cord2Rel[1],c,-s);
                    rotate(off[0],off[1],c,-s);
                    rotate(rot[0],rot[1],c,-s);
                    d -= cid - ped;
                    sourceIndex = cii;
                }

                //
                //d = 0;
                //...
            } else {//didn't hit a wall nor portal
                //printf("%i\n",iterations);
                vec3(E.cord1, d,0,0);
                E.cord1[0] -= off[0]; E.cord1[1] -= off[1];
                rotate(E.cord1[0],E.cord1[1],rot[0],-rot[1]);
                //
                vec3(E.cord2, E.radius*cord2Rel[0] + d, E.cord2[1] = E.radius*cord2Rel[1]*mirror,0);
                E.cord2[0] -= off[0]; E.cord2[1] -= off[1];
                rotate(E.cord2[0],E.cord2[1],rot[0],-rot[1]);
                //
                vec3(E.compass, d + dx,-dy*mirror,0);
                E.compass[0] -= off[0]; E.compass[1] -= off[1];
                rotate(E.compass[0],E.compass[1],rot[0],-rot[1]);
                //
                d = 0;
                //
                clearDuplicates(E);
                i=(int)paw[E.currentWorld];
                bool inWall = false;
                double cp[2];
                while(i < paw[E.currentWorld + 1] && !inWall){
                    type = ((char)paw[i+7])&7;
                    bool isWall = (type > 3);
                    if(!isWall || !E.bounce){
                        double p1[] = {paw[i],paw[i + 1]}, p2[] = {paw[i + 3],paw[i + 4]};
                        double p0[2] = {E.cord1[0],E.cord1[1]};
                        p2[0] -= p1[0]; p2[1] -= p1[1];
                        p0[0] -= p1[0]; p0[1] -= p1[1];
                        double r2 = p2[0]*p2[0] + p2[1]*p2[1];
                        r = sqrt(r2);
                        type &= 3;
                        double d2;
                        if(type == 0){
                            double al = paw[i + 6], ala = r2*al;//TODO remove all uses of al and just use ala
                            double c = 1 - al/2;
                            double s = sign(p0[1]*p2[0] - p0[0]*p2[1])*safe_sqrt(c*c - 1);
                            double endPt[2] = {p2[0],p2[1]};
                            rotate(endPt[0],endPt[1],c,s);
                            double C = r/sqrt(p0[0]*p0[0] + p0[1]*p0[1]);
                            cp[0] = C*p0[0]; cp[1] = C*p0[1];
                            d2 = (cp[0] - p2[0])*(cp[0] - p2[0]) + (cp[1] - p2[1])*(cp[1] - p2[1]);
                            if( d2 > ala && al != 4){
                                cp[0] = endPt[0]; cp[1] = endPt[1];
                            }
                        } else {
                            double C = (p0[0]*p2[0] + p0[1]*p2[1])/r2;
                            if(C < 0) C = 0;
                            if(C > 1) C = 1;
                            cp[0] = C*p2[0]; cp[1] = C*p2[1];
                        }
                        double tempv[2];
                        tempv[0] = cp[0] - p0[0]; tempv[1] = cp[1] - p0[1];
                        double disToCp = sqrt(tempv[0]*tempv[0] + tempv[1]*tempv[1]);
                        if( (disToCp < E.radius + wallThickness && isWall) || (disToCp < ped && !isWall)){
                            d = ped - disToCp;
                            if(isWall) d += wallThickness + E.radius;
                            else d += ped;
                            d*=1.01;
                            dx = E.compass[0] - E.cord1[0];
                            dy = E.compass[1] - E.cord1[1];
                            rotate(tempv[0],tempv[1],dx,-dy);
                            dx = -tempv[0]/disToCp; dy = -tempv[1]/disToCp;
                            entranceAng[0] = dx; entranceAng[1] = -dy;
                            type = -1;
                            sourceIndex = -1;
                            copypt(E.cord1,sourcept1);
                            copypt(E.compass,sourcept2);
                            inWall = true;
                            mirror = 1;
                        }
                        if(!inWall && disToCp < E.radius && !isWall){
                            double relLoc;
                            double sideSign = sign(p0[1]*p2[0] - p0[0]*p2[1]);
                            if(type == 0){
                                relLoc = r*sideSign*acos(1 - (d2/r2)/2);
                                disToCp *= sign(p0[0]*p0[0] + p0[1]*p0[1] - r2);
                            } else {
                                relLoc = sqrt(cp[0]*cp[0] + cp[1]*cp[1]) - r/2;
                                disToCp *= -sideSign;
                            }
                            double relPoses[3][3];
                            relPoses[0][0] = E.cord2[0] - E.cord1[0];
                            relPoses[0][1] = E.cord2[1] - E.cord1[1];
                            relPoses[1][0] = E.compass[0] - E.cord1[0];
                            relPoses[1][1] = E.compass[1] - E.cord1[1];
                            relPoses[2][0] = tempv[0];
                            relPoses[2][1] = tempv[1];
                            addDuplicate(E,i,relLoc,disToCp,relPoses);
                        }
                    }
                    i += sizeOfPow;
                }
            }
            //...
        } else if(iaunit > 0){//spherical
            double ds = sin(d*iaunit), dc = cos(d*iaunit);
            double cr = dot(sourcept1,sourcept2), sr = sqrt(1 - cr*cr);
            if(type < 4){
                    s2matto(sourcept1,transform);
                matxpt(transform,sourcept2);
                sourcept2[0] /= sr; sourcept2[1] /= sr;
                rotateXY(transform, sourcept2[0],-sourcept2[1]);
                if(sourceIndex > 0){
                    double a = -entranceLoc*iaunit/sr;
                    rotateXY(transform, cos(a),sin(a));
                    rotateXZ(transform, cr,sr);
                }
                rotateXY(transform, entranceAng[0],entranceAng[1]);
                if(sourceIndex == -1){
                    matxpt(transform,E.cord2,cord2Rel);
                    normalizeXY(cord2Rel);
                }
            }
            int i = (int)paw[E.currentWorld];
            double cip2[3],cip[3];
            double cidr = s2disrank(dc,ds);
            while(i<(int)paw[E.currentWorld + 1]){
                double p1[3]={paw[i],paw[i+1],paw[i+2]}, p2[3]={paw[i+3],paw[i+4],paw[i+5]};
                cr =dot(p1,p2);
                sr = 1 - cr*cr;//actually sr^2, just reusing space
                matxpt(transform,p1);
                matxpt(transform,p2);
                double y1y1=p1[1]*p1[1];
                if(y1y1 < sr){
                    double cipp[3], cipm[3];//closest intersection point plus and minus
                    double s=sqrt(sr-y1y1), D=1-y1y1;//angle limit, (cos(a)-1), decreases with angle
                    double al=paw[i+6], ala=1+sr*al, disrank;
                    cipp[0]=(p1[0]*cr+p1[2]*s)/D;
                    cipm[0]=(p1[0]*cr-p1[2]*s)/D;
                    cipp[2]=(p1[2]*cr-p1[0]*s)/D;
                    cipm[2]=(p1[2]*cr+p1[0]*s)/D;
                    cipp[1]=0;cipm[1]=0;

                    disrank = s2disrank(cipp[2],cipp[0]);
                    if( ( al == -2 || dot(cipp,p2) > ala ) && disrank < cidr && (sourceIndex != i || cipp[2] < cipm[2]) ){
                        cii=i;
                        copypt(p1,cip1);
                        copypt(p2,cip2);
                        copypt(cipp,cip);
                        cidr=disrank;
                    }
                    disrank=s2disrank(cipm[2],cipm[0]);
                    if( ( al == -2 || dot(cipm,p2) > ala ) && disrank < cidr && (sourceIndex !=i || cipm[2] < cipp[2]) ){
                        cii=i;
                        copypt(p1,cip1);
                        copypt(p2,cip2);
                        copypt(cipm,cip);
                        cidr=disrank;
                    }

                }
                i+=sizeOfPow;
            }//end of portal checks*/
            if(cii > 0){
                dat = int(paw[cii + 7]);
                type = dat&7;
                if(type > 3){
                    rotateXZ(transform,cip[2],cip[0]);
                    rotate(cip1[0],cip1[2],cip[2],cip[0]);
                    //
                    double r2 = cip1[0]*cip1[0] + cip1[1]*cip1[1];
                    double c = 2*cip1[1]*cip1[1]/r2 - 1, s = -2*cip1[0]*cip1[1]/r2;
                    //
                    rotate(cord2Rel[0],cord2Rel[1],c,-s);
                    rotateXY(transform,c,-s);
                    d -= acos(cip[2])/iaunit - ped;
                    sourceIndex = cii;
                } else {
                    E.currentWorld = ((dat)>>3)&511;
                    di = int(paw[E.currentWorld]) + sizeOfPow*((dat>>12)&1023);
                    sourceIndex = di;
                    sourcept1[0]=paw[di];sourcept1[1]=paw[di+1];sourcept1[2]=paw[di+2];
                    sourcept2[0]=paw[di+3];sourcept2[1]=paw[di+4];sourcept2[2]=paw[di+5];
                    if(((dat>>22)&1)==1) side = -1;
                    else side = 1;
                    if(((dat>>23)&1)==1) {mirror = -1; E.mirror *=-1;}
                    else mirror = 1;
                    //
                    double topbot = s2side(cip2,cip1,cip);
                    entranceLoc = dot(cip,cip2);
                    cr = dot(cip1,cip2);
                    sr = 1 - cr*cr;
                    entranceLoc = (entranceLoc - cr*cr)/sr;
                    if(entranceLoc > 1) entranceLoc = 1;
                    if(entranceLoc < -1) entranceLoc = -1;
                    sr = sqrt(sr);
                    entranceLoc = -acos(entranceLoc)*sr*mirror*topbot*side/iaunit;
                    rotate(cip1[0],cip1[2],cip[2],cip[0]);
                    entranceAng[0] = -cip1[0]/sr*side;
                    entranceAng[1] = -cip1[1]/sr*side*mirror;
                    d -= arctan(cip[0],cip[2])/iaunit;
                }
            } else {
                rotateXZ(transform, dc,ds);
                vec3(E.cord1, 0,0,1);
                double sinEr = sin(E.radius*iaunit);
                vec3(E.cord2, cord2Rel[0]*sinEr,mirror*cord2Rel[1]*sinEr,cos(E.radius*iaunit));
                vec3(E.compass, dx,-dy*mirror,0);
                double itransform[3][3];
                transpose(transform,itransform);
                matxpt(itransform,E.cord1);
                backOnSphere(E.cord1);
                matxpt(itransform,E.cord2);
                backOnSphere(E.cord2);
                matxpt(itransform,E.compass);
                backOnSphere(E.compass);
                d = 0;
                //
                clearDuplicates(E);
                double cosEr = cos(E.radius*iaunit);
                double cosErcol = cos( (E.radius + wallThickness)*iaunit);
                double cosped = cos(ped*iaunit);
                bool inWall = false;
                i = (int)paw[E.currentWorld];
                while(i < paw[E.currentWorld + 1] && !inWall){
                    bool isWall = (int(paw[i + 7])&4) != 0;
                    if(!isWall || !E.bounce){
                        double p0[3], p1[] = {paw[i],paw[i+1],paw[i+2]}, p2[]={paw[i+3],paw[i+4],paw[i+5]};
                        cr = dot(p1,p2);
                        sr = 1 - cr*cr;//actually sr^2
                        double limit = 1 + paw[i+6]*sr;
                        copypt(E.cord1,p0);
                        double d10 = dot(p1,p0);//You thought it would be a better name, but it was ME,...
                        double A = sqrt(sr/(1-d10*d10)), B = cr - A*d10;
                        double cp[3];
                        for(int n = 0; n< 3; n++) cp[n] = A*p0[n] + B*p1[n];
                        double relLoc = dot(cp,p2);
                        double topbot = s2side(p1,p2,cp);
                        if(relLoc < limit){
                            relLoc = limit;
                            int o;
                            if(topbot > 0) o = 8;
                            else o = 11;
                            for(int n = 0; n < 3; n++) cp[n] = paw[i + o + n];
                        }
                        double cosDisToCp = dot(p0,cp);
                        if( (isWall && cosDisToCp > cosErcol) || (!isWall && cosDisToCp > cosped)){
                            inWall = true;
                            s2matto(p0,transform);
                            double tv[3];
                            matxpt(transform,E.compass,tv);
                            matxpt(transform,cp);
                            rotate(cp[0],cp[1],-tv[0],tv[1]);
                            normalizeXY(cp);
                            dx = cp[0];
                            dy = cp[1];
                            entranceAng[0] = dx;
                            entranceAng[1] = -dy;
                            mirror = 1;
                            type = 0;
                            sourceIndex = -1;
                            copypt(E.cord1,sourcept1);
                            copypt(E.cord2,sourcept2);
                            d = ped - acos(cosDisToCp)/iaunit;
                            if(isWall) d += wallThickness + E.radius;
                            else d += ped;
                        }
                        if(!inWall && cosDisToCp > cosEr){
                            s2matto(p0,transform);
                            double relPoses[3][3];
                            matxpt(transform,E.cord2,relPoses[0]);
                            matxpt(transform,E.compass,relPoses[1]);
                            matxpt(transform,cp,relPoses[2]);
                            relLoc = acos( (relLoc - cr*cr)/sr )*sqrt(sr)/iaunit*topbot;
                            double disToCp = acos(cosDisToCp)/iaunit*sign(cr - d10);
                            addDuplicate(E,i,relLoc,disToCp,relPoses);
                        }
                    }
                    i += sizeOfPow;
                }
                //...
            }
            //...
        }
        /*
        if(cii > 0){
            E.currentWorld = ((dat)>>3)&511;
            di = int(paw[E.currentWorld]) + sizeOfPow*((dat>>12)&1023);
            vec3(sourcept1, paw[di],paw[di + 1],paw[di + 2]);
            vec3(sourcept2, paw[di + 3],paw[di + 4],paw[di + 5]);
            sourceIndex = di;
            if(((dat>>22)&1)==1) side=-1;
            else side=1;
            if(((dat>>23)&1)==1) {mirror = -1; E.mirror *=-1;}
            r *= side;
            entranceAng[0] = -cip1[0]/r; entranceAng[1] = -cip1[1]/r*mirror;
            entranceLoc *= mirror*side;
            type = char(int(paw[di+7])&3);
        }//*/
    }
    if(iterations == maxits) abortMove = true;
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

/*void renderFrame(int windowres, int windowScale, int height, int width, float zoom, unsigned int shaderProgram, unsigned int framestage1, unsigned int framebuffer1, unsigned int postprocessfull, GLFWwindow* window, double T){
    glfwMakeContextCurrent(window);
    printf("1\n");
    glDrawArrays(GL_TRIANGLES, 0, 3);
    printf("0\n");
            //
            //2nd pass
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glUseProgram(postprocessfull);
            glUniform1i(glGetUniformLocation(postprocessfull, "height"),height);
            glUniform1i(glGetUniformLocation(postprocessfull, "width"),width);
            glUniform1i(glGetUniformLocation(postprocessfull, "windowScale"),windowScale);
            glUniform1f(glGetUniformLocation(postprocessfull, "zoom"),zoom);
            glDrawArrays(GL_TRIANGLES, 0, 3);

            // check and call events and swap the buffers
            glfwSwapBuffers(window);
    glUseProgram(shaderProgram);
    //....
}//*/

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
        double p1[3] = {dc*sl,ds,0}, p2[3] = {sinh(1.0/9), 0, cosh(1.0/9)};
        p1[2] = sqrt(p1[0]*p1[0] + p1[1]*p1[1] + 1);
        double iso[3][3];
        h2matto(p1,iso);
        double isoi[3][3];
        h2invert(iso,isoi);
        matxpt(isoi,p2);
        std::vector<std::vector<float>> world0={
        {lc,0,sl, sl,0,lc, dc,mkdest(2,0,2,0,0)},
        {0,1.00,0, 0,ds,dc, limit,mkdest(2,0,3,0,0)},
        {lc,0,-sl, -sl,0,lc, dc,mkdest(2,0,0,0,0)},
        {0,1.00,0, 0,-ds,dc, limit,mkdest(2,0,1,0,0)},
        {p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], inf,4},
        {-p1[0],p1[1],p1[2], -p2[0],p2[1],p2[2], inf,4},
        {p1[0],-p1[1],p1[2], p2[0],-p2[1],p2[2], inf,4},
        {-p1[0],-p1[1],p1[2], -p2[0],-p2[1],p2[2], inf,4},
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
            {0,2,sqrt(5), 0,0,1, inf,mkdest(0,1,0,0,0)},//d*d/2 + 1 = L for horocycles
            {0,-3,sqrt(8), 0,-1,sqrt(2), 3,6}
        };
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(-1);
        //std::vector<std::vector<float>> world1={{2,3,0, -2,3,0, 4,mkdest(1,0,0,0,0)}};
        std::vector<std::vector<float>> world1={{2,0,0, 2,2,0, 4,mkdest(0,0,0,0,0)}};
        pawbuffer.push_back(world1);
        worldCurvatures.push_back(0);
        pw = 0; plp[0] = 1; duppw = -1;//signals that the player isn't in a portal
        vec3(pl, 0,0,1); vec3(camRef, 1,0,sqrt(2));
        vec3(pl2, prsh,0,prch);
        //
        //vec3(pl, 0,0,0);vec3(pl2, pl[0],pl[1]+pr,0);vec3(camRef, pl[0]+1,pl[1],0);pw=1;
    }//*/
    //NOTE: isr2 = 1/sqrt(2), sr2 = sqrt(2), I just got tired of typing it in testing

    /*{//test stage
        std::vector<std::vector<float>> world0={
            {0,2,0, 0.5,2,0, 2,mkdest(0,1,0,1,0)},//portal, type 0 (arc), linked to world1, index 0, other side = 1, mirror = 0
            {0,1.25,0, 1,1.25,0, 0,mkdest(1,0,1,0,1)},//portal, type 1(E2 line), linked to world 0, index 1 (itself), side = 0, mirror = 1
            {0,-5,0, 0,-4,0, 2,mkdest(0,0,3,0,0)},
            {3,-2,0, 2,-2,0, 2,mkdest(0,0,2,0,0)},
            {-2,5,0, -2,1,0, 0,mkdest(1,2,0,0,0)},
            {4,0,0, 4.5,0,0, 4,mkdest(0,3,0,0,0)},
            {-3,0,0, -2,0,0, 4,mkdest(0,4,0,0,0)},
            {0,4,0, sinh(0.5),4,0, 4,mkdest(0,2,2,0,0)},
            {1,1,0, 2,1,0, 0,5},//E2 line wall
            {2,2,0, 3,2,0, 2,4},//arc wall
            {-1.5,1,0, -2.5,1,0, 0,5},
            {0,1,0, 0,1.5,0, 0,5},
            {0,2.5,0, 1.0/9,2.5,0, 4,4},
            {1,1,0, 1,2,0, 0,5},
        };
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(0);//euclidian
        /////////////////////////
        std::vector<std::vector<float>> world1={
            {0,0,1, 0,0.5,0.866025403784, -1,mkdest(0,0,0,1,0)}, //portal, type 0 (arc), linked to world0, index 0, other side = 1, mirror = 0
            {0,0,1, 0,1,0, -0.25,4}, //arc wall
            {1,0,0, 0,0.923879533,-0.382683432, -0.076120467,4},
            {0,-1,0, -isr2,0,isr2, -0.1,4},
            {0.500000000,0.000000000,0.866025404, 0.592943927,0.000000000,0.805243751, -2,4},
        };
        pawbuffer.push_back(world1);
        worldCurvatures.push_back(1);//spherical
        /////////////////////
        //coments WIP
        std::vector<std::vector<float>> world2={
            {0,2,sqrt(4), 0,0,1, 3,mkdest(1,0,4,0,0)},//d*d/2 + 1 = -L for horocycles
            {0,1,1, 0,2,sqrt(5), 2,mkdest(1,2,1,0,1)},
            {1.500207901,3.017860025,3.515409347, 2.507112635,5.043377651,5.720251034, inf,mkdest(0,0,7,0,0)},
            {0,3,-sqrt(8), 0,-1,sqrt(2), 3,6},//hyperciclic wall
            {0.984237845,-1.983477923,-1.975578078, -2.004547132,2.004228044,3.005850910, 1.248541267,6},
        };
        pawbuffer.push_back(world2);
        worldCurvatures.push_back(-1);//hyperbolic
        //
        std::vector<std::vector<float>> world3={
         {0,0,-1, 0,0.5,0.866025403784, -2,mkdest(0,0,5,0,0)},
        };
        pawbuffer.push_back(world3);
        worldCurvatures.push_back(1);
        //
        std::vector<std::vector<float>> world4={
         {0,0,1, 0,1,0, -2,mkdest(0,0,6,0,0)},
        };
        pawbuffer.push_back(world4);
        worldCurvatures.push_back(1);
        ///
        pw = 0;//player in world 0
        vec3(pl, 0,0,0);//player at 0,0,0
        vec3(pl2, pl[0],pl[1]+pr,0);//placing the default facing of the player
        vec3(camRef, pl[0]+1,pl[1],0);
        plp[0] = 1;//setting this to negative 1 mirrors
        duppw = -1;//signals that the player isn't in a portal
    }//*/

    {//showcase
        std::vector<std::vector<float>> world0={
        {2,2,0, -2,2,0, 0,mkdest(1,0,1,0,0)},
        {2,-2,0, -2,-2,0, 0,mkdest(1,0,0,0,0)},
        {-3,0,0, -2,0,0, 1,mkdest(0,0,3,0,0)},
        {2,0,0, 3,0,0, 1,mkdest(0,0,2,0,0)},
        {-5,6,0, -3,6,0, 0,mkdest(1,0,5,0,0)},
        {-6,5,0, -6,3,0, 0,mkdest(1,0,4,0,0)},
        {-6,-5,0, -5,-6,0, 0,mkdest(1,1,0,1,0)},
        {-1,6,0, 1,6,0, 0,mkdest(1,0,8,0,1)},
        {-1,4,0, 1,4,0, 0,mkdest(1,0,7,0,1)},
        {-1,-4,0, 1,-4,0, 0,mkdest(1,0,9,0,1)},
        {8,8,0, 8,8.5,0, 2-2*cos(1.25),mkdest(0,18,1,0,0)},
        {-5,-.5,0, -5,.5,0, 0,mkdest(1,3,0,0,0)},
        };
        pawbuffer.push_back(world0);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world1={
        {-5,6,0, -6,5,0, 0,mkdest(1,0,6,1,0)},
        {-3,4,0, -1,4,0, 1,mkdest(0,1,2,1,0)},
        {3,4,0, 1,4,0, 1,mkdest(0,1,1,1,0)},
        {4,0,0, 4,2,0, 1,mkdest(0,1,4,1,0)},
        {4,0,0, 4,-2,0, 1,mkdest(0,1,3,1,0)},
        {-3,-5,0, -1,-5,0, 0,mkdest(1,1,5,1,0)},
        {2,-6,0, 2,-4,0, 1,mkdest(0,1,6,1,0)},
        {-6,2,0, -4,2,0, 1,mkdest(0,1,7,1,1)},
        {-5,-3,0, -5,-1,0, 0,mkdest(1,1,8,1,1)},
        {-5,-6,0, -6,-6,0, 0,mkdest(1,2,0,0,0)},
        {-6,6,0, -6,6.5,0, 2-2*cos(1.25),mkdest(0,18,2,0,0)},
        };
        pawbuffer.push_back(world1);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world2={
        {0.5,6,0, -0.5,6,0, 0,mkdest(1,1,9,0,0)},
        {9,0,0, 1,0,0, 0.5,mkdest(0,2,2,1,0)},
        {-9,0,0, -1,0,0, 0.5,mkdest(0,2,1,1,0)},
        {0.5,-6,0, -0.5,-6,0, 0,mkdest(1,4,0,0,0)},
        {0,7,0, 0,7.5,0, 2-2*cos(1.25),mkdest(0,18,3,0,0)},
        {0.5,-6,0, 0.6,-6,0, 4,4},
        {-0.5,-6,0, -0.6,-6,0, 4,4},
        };
        pawbuffer.push_back(world2);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world3={
        {0.5,6,0, -0.5,6,0, 0,mkdest(1,0,11,0,0)},
        {-1,1.5,0,0,1.5,0,0,mkdest(1,3,1,0,1)},
        {4,0,0, 4,1,0, 0,mkdest(1,21,1,0,0)},
        {-0.5,-1.25,0, -0.5,-1.75,0, 2-2*cos(1.25),mkdest(0,18,4,0,0)},
        {-1+.00,0,0, -2.00,0,0, 2.001,4},
        {-1,1-.00,0, -1,2.00,0, 0,5},
        {-1,0,0, -3,0,0, 2.001,4},
        {-1.0,-1,0, 0.0,-1,0, 0,5},
        {-1.00,-2,0, 0.00,-2,0, 0,5},
        {0,-1,0, 1,0,0, 0,5},
        {1,-1.00,0, 1,0.00,0, 0,5},
        {1,1,0,0,1,0,0,5},
        {0,1,0,0,2,0,0,5},
        {0,2,0,1,1,0,0,5},
        {0,-1,0,isr2,-1-isr2,0,(2-2*isr2)+0.001,4},
        {4,0,0, 4.1,0,0, 4,4},
        {4,1,0, 4,1.1,0, 4,4},
        };
        //-0.169333       1.418431        0.000000        0.830667        1.418431        0.000000        3.000000        356.000000      -0.830667       1.418431        0.000000        3.000000        356.000000      0.028987        -0.032000
        //visual position BUG
        pawbuffer.push_back(world3);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world4={
        {4,0,0, 4,1,0, 0,mkdest(1,2,3,0,0)},
        {0,-2,0, 0,-3,0, 4,mkdest(0,5,0,0,0)},
        {-4,0,0, -3,0,0, 2,mkdest(0,6,0,0,0)},
        {3,-2,0, 3+isr2,-2,0, 4,mkdest(0,7,0,0,0)},
        {-3,-2,0, -3+isr2,-2,0, 4,mkdest(0,8,0,0,0)},
        {0,4,0, 0,4-sin(pi/4),0, 2-2*cos(pi/4),mkdest(0,9,0,0,0)},
        {-8,-8,0, -8.5,-8,0, 2-2*cos(1.25),mkdest(0,18,5,0,0)},
        {3,-2,0, 3.5,-2,0, 2-2*cos(1.25),mkdest(0,18,6,0,0)},
        {0,-2,0, 0,-2.5,0, 2-2*cos(1.25),mkdest(0,18,7,0,0)},
        {-3,-2,0, -3,-2.5,0, 2-2*cos(1.25),mkdest(0,18,13,0,0)},
        {-0.5,3.5,0, 0.1 - 0.5,3.5,0, 4,4},
        {0.5,3.5,0, 0.1 + 0.5,3.5,0, 4,4},
        {4,0,0, 4,0.1,0, 4,4},
        {4,1,0, 4,1.1,0, 4,4},
        };
        pawbuffer.push_back(world4);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world5={
        {0,0,1, 0,1,0, -3,mkdest(0,4,1,0,0)},
        };
        pawbuffer.push_back(world5);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world6={
        {0,0,1, 1,0,0, -1,mkdest(0,4,2,0,0)},
        };
        pawbuffer.push_back(world6);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world7={
        {0,0,1, 0,isr2,isr2, -3,mkdest(0,4,3,0,0)},
        };
        pawbuffer.push_back(world7);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world8={
        {0,0,-1, 0,isr2,isr2, -3,mkdest(0,4,4,0,0)},
        };
        pawbuffer.push_back(world8);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world9={
        {0,0,1, 0,isr2,isr2, cos(pi/4)-1,mkdest(0,4,5,0,0)},
        {1,0,0, isr2,isr2,0, cos(pi/8)-1,mkdest(0,9,2,1,0)},
        {1,0,0, isr2,-isr2,0, cos(pi/8)-1,mkdest(0,9,1,1,0)},
        {1,0,0, isr2,0,-isr2, cos(pi/8)-1,mkdest(0,9,4,0,0)},
        {0,0,1, isr2,0,isr2, cos(pi/8)-1,mkdest(0,9,3,0,0)},
        {-0.577350269,-0.577350269,-0.577350269, -0.430918160,-0.430918160,-0.792855017, -3,mkdest(0,10,0,0,0)},
        //{-0.577350269,-0.577350269,-0.577350269, -0.672209398,-0.672209398,-0.310272541, -1.01,4},
        {-0.506819707,-0.831013078,0.229239331, -0.628517767,0.538215910,0.561506043, -0.019588431,4},
        {0.5,0.5,0.707106781, 0.547418791,0.547418791,0.632981307, -2,4},
        {-0.5,0.5,0.707106781, -0.547418791,0.547418791,0.632981307, -2,4},
        {-0.656685568,0.309852332,0.687572154, -0.844812768,0.398226150,0.357361525, -0.132,4},
        };
        pawbuffer.push_back(world9);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world10={
        {0.577350269,0.577350269,0.577350269, -0.430918160,-0.430918160,-0.792855017, -3,mkdest(0,9,5,0,0)},
        {1,0,0, 0,0,1, -0.02,mkdest(0,11,0,0,0)},
        };
        pawbuffer.push_back(world10);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world11={
        {1,0,0, 0,0,1, -0.02,mkdest(0,10,1,0,0)},
        {1,0,0, isr2,isr2,0, cos(pi/8)-1,mkdest(0,12,0,0,0)},
        };
        pawbuffer.push_back(world11);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world12={
        {1,0,0, isr2,-isr2,0, cos(pi/8)-1,mkdest(0,11,1,0,0)},
        {-1,0,0, -isr2,isr2,0, cos(pi/8)-1,mkdest(0,13,0,0,0)},
        {1,0,0, isr2,isr2,0, cos(7*pi/8)-1,4},
        {-1,0,0, -isr2,-isr2,0, cos(7*pi/8)-1,4},
        };
        pawbuffer.push_back(world12);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world13={
        {0,0,1, isr2,0,sqrt(1.5), 1+(1-cos(pi/8))*0.5, mkdest(0,12,1,0,0)},
        {0,1,1, 0,3.75/2,4.25/2, 1.5*1.5/2.0+1,mkdest(1,14,0,1,0)},
        {0,sinh(-2)*cosh(0.5)+cosh(2)*sinh(0.5),cosh(2)*cosh(0.5)+sinh(-2)*sinh(0.5), 0,sinh(-2),cosh(2), sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,8,0,0)},//manual lorentz boost
        {0,1,1, 0,1.5/2,2.5/2, 0.5*0.5/2.0+1,5},
        {0,-1,1, 0,-1.5/2,2.5/2, inf,5},
        };
        pawbuffer.push_back(world13);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world14={
        {-1.5,3,0, 1.5,3,0, 0,mkdest(1,13,1,1,0)},
        {-3,-3,0, -3,0,0, 4,mkdest(0,15,0,0,0)},
        {3,-3,0, 3,-2.5,0, 2*(1 - cos(1.25)),mkdest(0,16,2,0,0)},
        {-3,-3,0, -3.5,-3,0, 2-2*cos(1.25),mkdest(0,18,9,0,0)}
        };
        pawbuffer.push_back(world14);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world15={
        {0,0,1, 0,3,sqrt(10),inf,mkdest(0,14,1,0,0)},
        };
        pawbuffer.push_back(world15);
        worldCurvatures.push_back(-1);
        {
            double rc = cosh(3), rs = sinh(3), hc = cosh(2.42264), hs = sinh(2.42264);
            double x1 = sqrt((1 + rs*rs)/(1 + hs*hs)), z1 = sqrt(x1*x1 - 1);
            double l1 = (rc*x1 - sqrt(rs*rs - hs*hs)*z1);
            double l2 = (rc*hc - hs*hs)*1.01;
            double p1[3] = {0,-hs,hc}, p2[3] = {0,sinh(.5),cosh(0.5)};
            lorentz(p2[1],p2[2],hc,-hs);
            double p12[3], p22[3];
            copypt(p1,p12);
            copypt(p2,p22);
            lorentz(p1[0],p1[2],cosh(1),sinh(1));
            lorentz(p2[0],p2[2],cosh(1),sinh(1));
            std::vector<std::vector<float>> world16={
                {x1,0,z1, z1,0,x1, l1,0 + 1*mkdest(2,16,1,0,0)},
                {x1,0,-z1, -z1,0,x1, l1,0 + 1*mkdest(2,16,0,0,0)},
                {p1[0],-p1[1],p1[2], p2[0],-p2[1],p2[2], sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,14,2,0,0)},
                {p12[0],p1[1],p12[2], p22[0],p22[1],p22[2], sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,17,0,0,0)},
                {sinh(-2)*cosh(0.5)+cosh(2)*sinh(0.5),0,cosh(2)*cosh(0.5)+sinh(-2)*sinh(0.5), sinh(-2),0,cosh(2), sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,10,0,0)},//manual lorentz boost
                {0,-1,0, 0,hs,hc, l2,6},
                {0,-1,0, 0,-hs,hc, l2,6},
            };
            pawbuffer.push_back(world16);
            worldCurvatures.push_back(-1);
        }
        std::vector<std::vector<float>> world17={
            {0,0,1, 0,sin(0.5),cos(0.5), cos(1.25)-1,mkdest(0,16,3,0,0)},
            {1,0,0, cos(0.5),-sin(0.5),0, cos(1.25)-1,mkdest(0,18,0,0,0)},
        };
        pawbuffer.push_back(world17);
        worldCurvatures.push_back(1);
        {
            double x1 = sinh(3), z1 = cosh(3), x2 = x1, z2 = z1, A = 2*pi/32;
            int N = 0;
            lorentz(x2,z2,cosh(.5),-sinh(0.5));
            std::vector<std::vector<float>> world18={//
                {x1,0,z1, x2,0,z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,17,1,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,0,10,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,1,10,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,2,4,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,3,3,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,4,6,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,4,7,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,4,8,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,13,2,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,14,3,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,16,4,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,19,0,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,19,1,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,4,9,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,21,0,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,23,0,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,22,0,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,24,0,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,25,0,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,26,0,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,30,4,0,0)},
                {x1*cos(++N*A),x1*sin(N*A),z1, x2*cos(N*A),x2*sin(N*A),z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,33,2,0,0)},
            };
            pawbuffer.push_back(world18);
            worldCurvatures.push_back(-1);
        }
        std::vector<std::vector<float>> world19={
            {0,0,0, 0,-0.5,0, 2-2*cos(1.25),mkdest(0,18,11,0,0)},
            {0,6,0, 0,6.5,0, 2-2*cos(1.25),mkdest(0,18,12,0,0)},
            {-400,3,0, 400,3,0, 0,mkdest(1,20,0,0,0)}
        };
        pawbuffer.push_back(world19);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world20={
            {0,1,1, 0,0,1, inf,mkdest(1,19,2,0,0)}
        };
        pawbuffer.push_back(world20);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world21={
            {-6,5.5,0, -6,6,0, 2-2*cos(1.25),mkdest(0,18,14,0,0)},
            {3,1.5,0, 3,.5,0, 0,mkdest(1,3,2,0,0)},
            {3,1.5,0, 3,1.6,0, 4,4},
            {3,.5,0, 3,.6,0, 4,4},
            {0,0,0, 97,0,0, 0,5},
            {97,0,0, 97,-97,0, 0,5},
            {0,-97,0, 97,-97,0, 0,5},
            {0,-97,0, 0,-63,0, 0,5},
            {-4,-63,0, 0,-63,0, 0,5},
            {-4,-63,0, -4,-85,0, 0,5},
            {-12.5,-85,0, -4,-85,0, 0,5},
            {-12.5,-85,0, -12.5,-91.5,0, 0,5},
            {-6.26,-91.5,0, -12.5,-91.5,0, 0,5},
            {-6.26,-91.5,0, -6.25,-103.5,0, 0,5},
            {103.5,-103.5,0, -6.25,-103.5,0, 0,5},
            {103.5,-103.5,0, 103.5,188.5,0, 0,5},
            {97,188.5,0, 103.5,188.5,0, 0,5},
            {97,188.5,0, 97,4,0, 0,5},
            {-4,4,0, 97,4,0, 0,5},
            {-4,4,0, -4,-31,0, 0,5},
            {0,-31,0, -4,-31,0, 0,5},
            {0,-31,0, 0,0,0, 0,5},
            {97,219.75,0, 103.5,219.75,0, 0,5},
            {103.5,232.75,0, 103.5,219.75,0, 0,5},
            {103.5,232.75,0, -235,231.75,0, 0,5},
            {-235.5,-101.25,0, -235,231.75,0, 0,5},
            {-235.5,-101.25,0, -63,-101.25,0, 0,5},
            {-63,-91.5,0, -63,-101.25,0, 0,5},
            {-63,-91.5,0, -43.75,-91.5,0, 0,5},
            {-43.75,-85,0, -43.75,-91.5,0, 0,5},
            {-43.75,-85,0, -66,-85,0, 0,5},
            {-66,-80.5,0, -66,-85,0, 0,5},
            {-66,-80.5,0, -70,-80.5,0, 0,5},
            {-70,-94.75,0, -70,-80.5,0, 0,5},
            {-70,-94.75,0, -93,-94.75,0, 0,5},
            {-93,-27.75,0, -93,-94.75,0, 0,5},
            {-93,-27.75,0, -70,-27.75,0, 0,5},
            {-70,-49.25,0, -70,-27.75,0, 0,5},
            {-70,-49.25,0, -66,-49.25,0, 0,5},
            {-66,-23.75,0, -66,-49.25,0, 0,5},
            {-66,-23.75,0, -101.75,-23.75,0, 0,5},
            {-101.75,-27.75,0, -101.75,-23.75,0, 0,5},
            {-101.75,-27.75,0, -97,-27.75,0, 0,5},
            {-97,-94.75,0, -97,-27.75,0, 0,5},
            {-229,-94.75,0, -97,-94.75,0, 0,5},
            {-229,-94.75,0, -229,-27.75,0, 0,5},
            {-133,-27.75,0, -229,-27.75,0, 0,5},
            {-133,-27.75,0, -133,18.25,0, 0,5},
            {-170,18.25,0, -133,18.25,0, 0,5},
            {-170,18.25,0, -170,14.25,0, 0,5},
            {-141.75,14.25,0, -170,14.25,0, 0,5},
            {-141.75,14.25,0, -141.75,14.25,0, 0,5},
            {-141.75,-23.75,0, -141.75,14.25,0, 0,5},
            {-229,-23.75,0, -141.75,-23.75,0, 0,5},
            {-229,14.25,0, -229,-23.75,0, 0,5},
            {-207,14.25,0, -229,14.25,0, 0,5},
            {-207,18.25,0, -207,14.25,0, 0,5},
            {-229,18.25,0, -207,18.25,0, 0,5},
            {-229,226.25,0, -229,18.25,0, 0,5},
            {-50,226,0, -229,226.25,0, 0,5},
            {-50,63.16,0, -50,226,0, 0,5},
            {-94.91,18.25,0, -50,63.16,0, 0,5},
            {-100.75,18.25,0, -94.91,18.25,0, 0,5},
            {-100.75,14.25,0, -100.75,18.25,0, 0,5},
            {-93.25,14.25,0, -100.75,14.25,0, 0,5},
            {-46,61.5,0, -93.25,14.25,0, 0,5},
            {-46,226.25,0, -46,61.5,0, 0,5},
            {97,226.25,0, -46,226.25,0, 0,5},
            {97,219.75,0, 97,226.25,0, 0,5},
        };
        for(int n = 4; n < world21.size(); n++){
            for(int k = 0; k < 6; k++){
                world21[n][k] /= 81.0;
            }
            world21[n][1]--;
            world21[n][4]--;
        }
        pawbuffer.push_back(world21);
        worldCurvatures.push_back(0);
        std::vector<std::vector<float>> world22 = {
            {0,-sinh(2.5),cosh(2.5), 0,-sinh(3),cosh(3), sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,16,0,0)},
        };
        for(int n = 0; n < world21.size(); n++){
            if((int(world21[n][7])&7) == 5){
                double m1[3][3], m2[3][3];
                double p3[] = {world21[n][0]+.5,world21[n][1]+.5,0}, p4[] = {world21[n][3]+.5,world21[n][4]+.5,0}, p2[3] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                double d,s,c;
                double* p;
                p = &p2[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                p[2] = d*d/2 + 1;
                s = d/2*sqrt(d*d + 4)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p3[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                p[2] = d*d/2 + 1;
                s = d/2*sqrt(d*d + 4)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p4[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                p[2] = d*d/2 + 1;
                s = d/2*sqrt(d*d + 4)/d;
                p[0] *= s;
                p[1] *= s;
                double p2c[3];
                copypt(p2,p2c);
                std::vector<float> wall = h2threepointstowall(p2c,p3,p4);
                world22.push_back(wall);
            }
        }
        pawbuffer.push_back(world22);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world23 = {
            {0,-sinh(2.5),cosh(2.5), 0,-sinh(3),cosh(3), sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,15,0,0)},
        };
        for(int n = 0; n < world21.size(); n++){
            if((int(world21[n][7])&7) == 5){
                double m1[3][3], m2[3][3];
                double p3[] = {world21[n][0]+.5,world21[n][1]+.5,0}, p4[] = {world21[n][3]+.5,world21[n][4]+.5,0}, p2[3] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                double d,s,c;
                double* p;
                p = &p2[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                p[2] = cosh(d);
                s = sinh(d)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p3[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                p[2] = cosh(d);
                s = sinh(d)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p4[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                p[2] = cosh(d);
                s = sinh(d)/d;
                p[0] *= s;
                p[1] *= s;
                double p2c[3];
                copypt(p2,p2c);
                std::vector<float> wall = h2threepointstowall(p2c,p3,p4);
                world23.push_back(wall);
            }
        }
        pawbuffer.push_back(world23);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world24  = {
            {0,-sinh(2.5),cosh(2.5), 0,-sinh(3),cosh(3), sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,17,0,0)},
        };
        for(int n = 0; n < world21.size(); n++){
            double scale = 3;
            if((int(world21[n][7])&7) == 5){
                double m1[3][3], m2[3][3];
                double p3[] = {world21[n][0]+0.75,world21[n][1]+.2,0}, p4[] = {world21[n][3]+0.75,world21[n][4]+.2,0}, p2[3] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                double d,s,c;
                double* p;
                p = &p2[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = atanh(d/scale);
                p[2] = cosh(s);
                s = sinh(s)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p3[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = atanh(d/scale);
                p[2] = cosh(s);
                s = sinh(s)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p4[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = atanh(d/scale);
                p[2] = cosh(s);
                s = sinh(s)/d;
                p[0] *= s;
                p[1] *= s;
                double p2c[3];
                copypt(p2,p2c);
                std::vector<float> wall = h2threepointstowall(p2c,p3,p4);
                world24.push_back(wall);
            }
        }
        pawbuffer.push_back(world24);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world25 = {
                {0,-sinh(2.5),cosh(2.5), 0,-sinh(3),cosh(3), sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,18,0,0)},
        };
        for(int n = 0; n < world21.size(); n++){
            double scale = 3;
            if((int(world21[n][7])&7) == 5){
                double m1[3][3], m2[3][3];
                double p3[] = {world21[n][0]+0.75,world21[n][1]+.2,0}, p4[] = {world21[n][3]+0.75,world21[n][4]+.2,0}, p2[3] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                double d,s,c;
                double* p;
                p = &p2[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = 2*atanh(d/scale);
                p[2] = cosh(s);
                s = sinh(s)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p3[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = 2*atanh(d/scale);
                p[2] = cosh(s);
                s = sinh(s)/d;
                p[0] *= s;
                p[1] *= s;
                p = &p4[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = 2*atanh(d/scale);
                p[2] = cosh(s);
                s = sinh(s)/d;
                p[0] *= s;
                p[1] *= s;
                double p2c[3];
                copypt(p2,p2c);
                std::vector<float> wall = h2threepointstowall(p2c,p3,p4);
                world25.push_back(wall);
            }
        }
        pawbuffer.push_back(world25);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world26;
        {
            double y1 = -sinh(2.5), z1 = cosh(2.5), y2 = y1, z2 = z1;
            lorentz(y2,z2,cosh(.5),-sinh(.5));
            world26 = {
                {0,y1,z1, 0,y2,z2, sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,19,0,0)},
            };
        }
        for(int n = 0; n < world21.size(); n++){
            if((int(world21[n][7])&7) == 5){
                double m1[3][3], m2[3][3];
                double p3[] = {world21[n][0]+.35,world21[n][1]+.5,0}, p4[] = {world21[n][3]+.35,world21[n][4]+.5,0}, p2[3] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                double d,s,c,a;
                double* p;
                p = &p2[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = sinh(d);
                c = cosh(d);
                a = arctan(p[0],p[1]);
                if(a > pi) a -= 2*pi;
                a = d*a/s;
                p[0] = s*sin(a);
                p[1] = s*cos(a);
                p[2] = c;
                p = &p3[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = sinh(d);
                c = cosh(d);
                a = arctan(p[0],p[1]);
                if(a > pi) a -= 2*pi;
                a = d*a/s;
                p[0] = s*sin(a);
                p[1] = s*cos(a);
                p[2] = c;
                p = &p4[0];
                d = sqrt(p[0]*p[0] + p[1]*p[1]);
                s = sinh(d);
                c = cosh(d);
                a = arctan(p[0],p[1]);
                if(a > pi) a -= 2*pi;
                a = d*a/s;
                p[0] = s*sin(a);
                p[1] = s*cos(a);
                p[2] = c;
                double p2c[3];
                copypt(p2,p2c);
                std::vector<float> wall = h2threepointstowall(p2c,p3,p4);
                world26.push_back(wall);
            }
        }
        pawbuffer.push_back(world26);
        worldCurvatures.push_back(-1);
        std::vector<std::vector<float>> world27;
        {
            for(int n = 0; n < 7*0 + 1*world21.size(); n ++){
                if((int(world21[n][7])&7) == 5){
                    double p3[] = {world21[n][0]+0.75,world21[n][1]+.2,0}, p4[] = {world21[n][3]+0.75,world21[n][4]+.2,0}, p5[] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                    double *p;
                    double d;
                    p = &p3[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(d);
                    d = sin(d)/d;
                    p[0] *= d; p[1] *= d;
                    p = &p4[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(d);
                    d = sin(d)/d;
                    p[0] *= d; p[1] *= d;
                    p = &p5[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(d);
                    d = sin(d)/d;
                    p[0] *= d; p[1] *= d;
                    std::vector<float> wall = s2threepointstowall(p5,p3,p4);
                    world27.push_back(wall);
                }
            }
        }
        pawbuffer.push_back(world27);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world28;
        {
            double scale = 1;
            for(int n = 0; n < 7*0 + 1*world21.size(); n ++){
                if((int(world21[n][7])&7) == 5){
                    double p3[] = {world21[n][0]+0.75,world21[n][1]+1,0}, p4[] = {world21[n][3]+0.75,world21[n][4]+1,0}, p5[] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                    double *p;
                    double d;
                    p = &p3[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(atan(scale*d));
                    d = sin(atan(scale*d))/d;
                    p[0] *= d; p[1] *= d;
                    p = &p4[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(atan(scale*d));
                    d = sin(atan(scale*d))/d;
                    p[0] *= d; p[1] *= d;
                    p = &p5[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(atan(scale*d));
                    d = sin(atan(scale*d))/d;
                    p[0] *= d; p[1] *= d;
                    std::vector<float> wall = s2threepointstowall(p5,p3,p4);
                    world28.push_back(wall);
                }
            }
        }
        pawbuffer.push_back(world28);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world29;
        {
            double scale = 1;
            for(int n = 0; n < 7*0 + 1*world21.size(); n ++){
                if((int(world21[n][7])&7) == 5){
                    double p3[] = {world21[n][0]+0.4,world21[n][1]+.9,0}, p4[] = {world21[n][3]+0.4,world21[n][4]+.9,0}, p5[] = {(p3[0]+p4[0])/2,(p3[1]+p4[1])/2,0};
                    double *p;
                    double d;
                    p = &p3[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(2*atan(d/scale));
                    d = sin(2*atan(d/scale))/d;
                    p[0] *= d; p[1] *= d;
                    p = &p4[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(2*atan(d/scale));
                    d = sin(2*atan(d/scale))/d;
                    p[0] *= d; p[1] *= d;
                    p = &p5[0];
                    d = sqrt(p[0]*p[0] + p[1]*p[1]);
                    p[2] = cos(2*atan(d/scale));
                    d = sin(2*atan(d/scale))/d;
                    p[0] *= d; p[1] *= d;
                    std::vector<float> wall = s2threepointstowall(p5,p3,p4);
                    world29.push_back(wall);
                }
            }
        }
        pawbuffer.push_back(world29);
        worldCurvatures.push_back(1);
        std::vector<std::vector<float>> world30;
        {//psudosphere case c, with fun added
            double d = 1.5, ds = sinh(d), dc = cosh(d), l = pi/ds, sl = sinh(l), lc = cosh(l);
            double limit = dc*dc*lc - ds*ds;
            double p1[3] = {dc*sl,ds,0}, p2[3] = {sinh(1.0/9), 0, cosh(1.0/9)};
            p1[2] = sqrt(p1[0]*p1[0] + p1[1]*p1[1] + 1);
            double iso[3][3];
            h2matto(p1,iso);
            double isoi[3][3];
            h2invert(iso,isoi);
            matxpt(isoi,p2);
            std::vector<std::vector<float>> world30={
                {lc,0,sl, sl,0,lc, dc,mkdest(2,30,2,0,0)},
                {0,1.00,0, 0,ds,dc, limit,mkdest(2,30,3,0,1)},
                {lc,0,-sl, -sl,0,lc, dc,mkdest(2,30,0,0,0)},
                {0,1.00,0, 0,-ds,dc, limit,mkdest(2,30,1,0,1)},
                {0,sinh(.7),cosh(.7), 0,sinh(0.2),cosh(0.2), sinh(.5)*sinh(.5)*(1 - cos(1.25)) + 1,mkdest(0,18,20,0,0)},
                {p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], inf,4},
                {-p1[0],p1[1],p1[2], -p2[0],p2[1],p2[2], inf,4},
                {p1[0],-p1[1],p1[2], p2[0],-p2[1],p2[2], inf,4},
                {-p1[0],-p1[1],p1[2], -p2[0],-p2[1],p2[2], inf,4},
            };
            pawbuffer.push_back(world30);
            worldCurvatures.push_back(-1);
        }
        std::vector<std::vector<float>> world31 = {
            {1,0,1, sinh(2),0,cosh(2), 4,mkdest(1,31,1,1,0)},
            {1,0,1, -sinh(2),0,cosh(2), 4,mkdest(1,31,0,1,0)},
        };
        pawbuffer.push_back(world31);
        worldCurvatures.push_back(-1);
        {
            double scale = 2, rad = pr/scale;
            double s = sinh(pi), c = cosh(pi);
            double p1[3] = {0,1,sr2}, p2[3];
            copypt(p1,p2);
            lorentz(p2[1],p2[2],cosh(rad),-sinh(rad));
            lorentz(p1[0],p1[2],c,s);
            lorentz(p2[0],p2[2],c,s);
            //
            std::vector<std::vector<float>> world32 = {
                {0,1,0, 0,1,sr2, 2*c-1,mkdest(2,33,0,1,0)},
                {0,1,0, 0,-1,sr2, 2*c-1,mkdest(2,33,1,0,0)},
                {c,0,s, s,0,c, cosh(acosh(sr2)-rad),mkdest(2,32,3,1,0)},
                {-c,0,s, -s,0,c, cosh(acosh(sr2)-rad),mkdest(2,32,2,1,0)},
                //{p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], inf,4},
                {p1[0],-p1[1],p1[2], p2[0],-p2[1],p2[2], inf,4},
                {-p1[0],-p1[1],p1[2], -p2[0],-p2[1],p2[2], inf,4},
                {-p1[0],p1[1],p1[2], -p2[0],p2[1],p2[2], inf,4},
            };
            pawbuffer.push_back(world32);
            worldCurvatures.push_back(-1/scale);
            std::vector<std::vector<float>> world33 = {
                {1+scale*sr2,0,0, 1,0,0, 4, mkdest(0,32,0,1,0)},
                {-1-scale*sr2,0,0, -1,0,0, 4, mkdest(0,32,1,0,0)},
                {0,2,0, 0,2.5,0, 2-2*cos(1.25),mkdest(0,18,21,0,0)},
                {1+2*scale*sr2,0,0, 1+2*scale*sr2+pr,0,0, 4,4},
                {-1-2*scale*sr2,0,0, -1-2*scale*sr2+pr,0,0, 4,4},
            };
            pawbuffer.push_back(world33);
            worldCurvatures.push_back(0);
        }
        pw = 33;//player world
        vec3(pl,0,0,0);//player location
        //vec3(pl, -1/sqrt(3),-1/sqrt(3),-1/sqrt(3));
        vec3(pl2, pl[0],pl[1]+pr,0);//placing the default facing of the player
        vec3(camRef, pl[0]+1,pl[1],0);
        //vec3(camRef, 1,0,sr2);
        plp[0] = 1;//setting this to negative 1 mirrors
        duppw = -1;//signals that the player isn't in a portal
    }//*/
    {//3.0.0 testing

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
            if(pawbuffer[w][i].size() < 14) addExtraToPOW(pawbuffer[w][i], worldCurvatures[w]);
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

//considering chaning to
//msiiiiiiiwwwwwwwwwwwwttt
//321098765432109876543210

    int defaultController;
    int defaultWidth;
    int defaultHeight;
    int windowScale;
    double framesPerTick;
    {
        std::ifstream f("config.txt");
        if (!f.is_open()) {
            printf("ERROR:config.txt wasn't able to be openned\n");
            return 1;
        }
        std::string s;
        std::getline(f, s);
        defaultController = stoi(s);
        std::getline(f, s);
        defaultWidth = stoi(s);
        std::getline(f, s);
        defaultHeight = stoi(s);
        std::getline(f, s);
        windowScale = stoi(s);
        std::getline(f, s);
        framesPerTick = stoi(s)/tickRate;
        f.close();
    }
    //printf("%lf\n",framesPerTick);


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

    int width, height, windowres = std::min(defaultWidth,defaultHeight);;

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    GLFWwindow* window = glfwCreateWindow(defaultWidth, defaultHeight, "NightGun", NULL, NULL);
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
    glViewport(0, 0, defaultWidth, defaultHeight);
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
    float deadzone = 0.1;
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
    float zoom = 1;
    int64_t tickCount = 0;
    double framesToRender = 1;
    float dx = 0,dy = 0;
    float facingAngle[5]={0,0,0,1,0};//3rd one is just working space
    GLsync fence;
    GLint frameStatus = GL_SIGNALED;


    entity player;
    player.currentWorld = pw;
    copypt(pl,player.cord1);
    copypt(pl2,player.cord2);
    copypt(camRef,player.compass);
    player.mirror = plp[0];
    player.radius = pr;
    player.bounce = false;
    player.intersectingPortalIndex = -1;
    {
        entity playerDupe;
        playerDupe.intersectingPortalIndex = -1;
        playerDupe.currentWorld = -1;
        player.duplicate = &playerDupe;
    }

    int tickCountLooped = 0;
    double framesRendered = 0;

    entity entityLog[120], entityLog2[120];
    double dxLog[120], dyLog[120];
    int logIndex = 0;
    bool kHeld = false, lHeld = false;
    bool worldchange = false;

    while(!glfwWindowShouldClose(window)){
        //double deltaTime=glfwGetTime();
        double time = glfwGetTime();
        double NextTickTime = time+1/tickRate;

        // input
        //x:-0.0024175168946385384        dy:-0.0000000000000000000       pos:-0.7051436801807300370, 0.7006944950205539202, 0.1086260325569532031        cam:0.6988327265491615092, 0.6608337933806174291, 0.2737362194340526300         cw:12
        if(tickCount>0){
            double playerspeed = ps;
            bool running = false;
            int kx = 0, ky = 0;
            if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, true);
            if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) ky+=1;
            if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) kx-=1;
            if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) ky-=1;
            if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) kx+=1;//*/
            if(glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
                printE(player);
            }
            int logIndexp = logIndex;
            if(glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) {
                if(!kHeld) {logIndex--;if(logIndex < 0) logIndex = 120-1;}
                kHeld = true;
            } else kHeld = false;
            if(glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
                if(!lHeld) {logIndex++;if(logIndex >= 120) logIndex = 0;}
                lHeld = true;
            } else lHeld = false;
            if(glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS) {
                lHeld = false;kHeld = false;
            }
            if(logIndexp != logIndex){
                printf("\n%i\n",logIndex);
                copyEntityState(entityLog[logIndex],player);
                copyEntityState(entityLog2[logIndex],player.duplicate[0]);
                printE(entityLog[logIndex]);
                printE(entityLog2[logIndex]);
                printf("%lf\t%lf\n",dxLog[logIndex],dyLog[logIndex]);
            }
            if(glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) player.radius = pr;
            if(glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) player.radius = pr/10;
            player.duplicate[0].radius = player.radius;
            if(glfwGetKey(window, GLFW_KEY_PAGE_UP) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_PAGE_DOWN) == GLFW_PRESS){
                if(!worldchange){
                    if(glfwGetKey(window, GLFW_KEY_PAGE_UP) == GLFW_PRESS ) player.currentWorld++;
                    else player.currentWorld--;
                    double iaunit = worldCurvatures[player.currentWorld];
                    player.cord1[0] = 0; player.cord1[1] = 0;
                    if(iaunit != 0) player.cord1[2] = 1;
                    else player.cord1[2] = 0;
                    player.compass[0] = 1; player.compass[1] = 0;
                    if(iaunit < 0){
                        player.compass[2] = sr2;
                        h2dadtopt(-player.radius*iaunit, 1,0, player.cord2);
                    } else if(iaunit > 0){
                        player.compass[2] = 0;
                        s2dadtopt(player.radius*iaunit, 1,0, player.cord2);
                    } else {
                        player.compass[2] = 0;
                        vec3(player.cord2, player.radius,0,0);
                    }

                }
                worldchange = true;
            } else worldchange = false;

            /*if(glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {//with test stage
                vec3(pl,0.4356526867844326967,  0.2435568025340280884,  0.8665372414041405680);
                vec3(camRef,0.8999781773730995971,  -0.1346341949234345636, -0.4146237030053032124);
                pw = 1;
                iopip = 202;
                vec3(playerduploc,0.2445038746919800110,  2.4373467901204493025,   0);
                duppw = 0;
                diopip = 6;
                moveEntity(0.0170011930167675018,-0.0266319829970598221,pl,pl2,camRef,&pw,plp);
            }*/


            facingAngle[2] = 0;
            int controllerNum = GLFW_JOYSTICK_1+defaultController;
            for(int jid = GLFW_JOYSTICK_1; jid <= GLFW_JOYSTICK_LAST; jid++){
                if(glfwJoystickIsGamepad(jid)){
                    controllerNum = jid;
                    break;
                }
            }
            const float *axes = glfwGetJoystickAxes(controllerNum, &axesCount);
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
                const unsigned char* buttons = glfwGetJoystickButtons(controllerNum, &buttonCount);
                for(int n = 0; n < buttonCount; n++){
                    //if(buttons[n] == GLFW_PRESS) printf("BUTTON: %i\n",n);
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
                    facingAngle[4]=facingAngle[1]/facingAngle[2]*player.mirror;
                    //printf("%lf\t%lf\t%lf\n",facingAngle[2],facingAngle[3],facingAngle[4]);
                    float iaunit = worldCurvatures[player.currentWorld];
                    if(iaunit == 0){//to do: do something similar for the duplicate
                        double crc[2]={player.compass[0]-player.cord1[0],player.compass[1]-player.cord1[1]};
                        player.cord2[0]=pr*facingAngle[3];player.cord2[1]=pr*facingAngle[4];
                        rotate(player.cord2[0],player.cord2[1],crc[0],crc[1]);
                        player.cord2[0]+=player.cord1[0];player.cord2[1]+=player.cord1[1];
                        player.cord2[2]=0;
                    } else if(iaunit > 0){
                        double rot[3][3],roti[3][3];
                        s2matto(player.cord1,rot);
                        matxpt(rot,player.compass,player.cord2);//reusing memory space that will be redefined anyway
                        rotateXY(rot,player.cord2[0],-player.cord2[1]);
                        s2dadtopt(pr*iaunit,facingAngle[3],facingAngle[4],player.cord2);
                        transpose(rot,roti);
                        matxpt(roti,player.cord2);
                        backOnSphere(player.cord2);
                    } else if(iaunit < 0){
                        double rot[3][3],roti[3][3];
                        h2matto(player.cord1,rot);
                        matxpt(rot,player.compass,player.cord2);//reusing memory space that will be redefined anyway
                        rotateXY(rot,player.cord2[0],-player.cord2[1]);
                        h2dadtopt(-pr*iaunit,facingAngle[3],facingAngle[4],player.cord2);
                        h2invert(rot,roti);
                        matxpt(roti,player.cord2);
                        backOnHyperboloid(player.cord2);
                    }
                }
            }

            bool playerMoved = true;
            if(dy!=0||dx!=0) {
                if(running){
                    dx *= 2;
                    dy *= 2;
                }
                //pl,cref, pw,iopip, playerduploc,duppw,diopip
                logIndex++;if(logIndex >= 120) logIndex = 0;
                copyEntityState(player,entityLog[logIndex]);
                copyEntityState(player.duplicate[0],entityLog2[logIndex]);
                dxLog[logIndex] = dx;
                dyLog[logIndex] = dy;
                //...
                //MOVE
                entity playerBU;
                copyEntityState(player,playerBU);
                entity playerDupBU;
                copyEntityState(player.duplicate[0],playerDupBU);
                //

                bool abortMove = false;
                double displacement = sqrt(dx*dx + dy*dy);
                dx /= displacement; dy /= displacement;
                double maxStep = player.radius/2;
                double dxM = dx*maxStep, dyM = dy*maxStep;
                while(displacement > ped){
                    entity pbu;
                    entity pdbu;
                    if(!abortMove){
                        copyEntityState(player,pbu);
                        copyEntityState(player.duplicate[0],pdbu);
                    }
                    abortMove = false;
                    if(displacement > maxStep) moveEntity(dxM,dyM,player,abortMove);
                    else {
                        if(displacement < ped) displacement = ped;
                        moveEntity(dx*displacement,dy*displacement,player,abortMove);
                    }
                    if(abortMove){//if there was an error, revert and try again going half as far
                        copyEntityState(pbu,player);
                        copyEntityState(pdbu,player.duplicate[0]);
                        dxM /= 2;
                        dxM /= 2;
                        maxStep /= 2;
                        displacement /= 2;
                    } else displacement -= maxStep;
                }
                /*
                int steps = ceil(displacement/pr/2);
                dx /= steps; dy /= steps;
                while(steps > 0){
                    moveEntity(dx,dy,player,abortMove);

                    if(abortMove) steps = 0;
                    else steps--;
                }//*/
                if( (!player.bounce && player.currentWorld == playerBU.currentWorld && disSquaredAprox(player.cord1,playerBU.cord1,worldCurvatures[player.currentWorld]) < 0.000002)
                    || (playerBU.intersectingPortalIndex != -1 && player.currentWorld == playerDupBU.currentWorld && disSquaredAprox(player.cord1,playerDupBU.cord1,worldCurvatures[player.currentWorld]) < 0.000002)
                    ){//Without this, wedging yourself in a corner in non euclidian geometry rotates the camera due to holonomy
                    //printf("%lf\n",disSquaredAprox(player.cord1,playerBU.cord1,worldCurvatures[player.currentWorld]));
                    copyEntityState(playerBU,player);
                    copyEntityState(playerDupBU,player.duplicate[0]);

                    playerMoved = false;
                    //printf("aborted\n");
                }//*/
            } else playerMoved = false;

            if(!playerMoved && player.intersectingPortalIndex >= 0 && facingAngle[2] > 0){
                updateDuplicateRot(player,facingAngle[3],facingAngle[4]);
                //updateDuplicateRot(pl,pl2,playerduploc,playerduploc2,pw,iopip);
            }
            facingAngle[2] = 0;

            if(glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && zoom < 1) zoom+=pow(2,floor(log2(zoom)-4));
            if(glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS && zoom>pr/pi) zoom-=pow(2,floor(log2(zoom)-4));
        }



        //Draw
        int widthp=width, heightp=height;
        glfwGetFramebufferSize(window, &width, &height);
        if(width!=widthp||height!=heightp) {
            windowres = std::min(width,height);
            printf("window dim:\n%i\tx\t%i\n",width,height);
            /*for(int n=0;n<paw[(int)paw[0]-1];n++) {
                if((int)paw[n]==paw[n]) printf("%i\t%i\n",n,(int)paw[n]);
                else printf("%i\t%f\n",n,paw[n]);
            }//*/
        }

        float T=1+sin(17.2787595947*glfwGetTime())/2;

        // rendering commands here

        framesToRender += framesPerTick;

        //double deltaTime=glfwGetTime();

        glGetSynciv(fence, GL_SYNC_STATUS, sizeof(GLint), NULL, &frameStatus);
        //if((tickCount&63) == 0) printf("%i\n",frameStatus);
        if(framesToRender >= 1 && frameStatus == GL_SIGNALED){
            framesRendered++;
            glfwSwapBuffers(window);
            framesToRender -= int(framesToRender);
            {
                for(int w=0;w<worldCurvatures.size();w++){
                    worldCurvs[w]=worldCurvatures[w];
                }
                glUniform1fv(glGetUniformLocation(shaderProgram, "iaunits"),128,worldCurvs);
                shader_data[0]=pawbuffer.size()+1;//serializing
                for(int w=0;w<pawbuffer.size();w++){
                    shader_data[w+1]=shader_data[w]+pawbuffer[w].size()*sizeOfPow;
                    for(int i=0;i<pawbuffer[w].size();i++){
                        for(int p=0;p<sizeOfPow;p++){
                            shader_data[(int)shader_data[w]+sizeOfPow*i+p]=pawbuffer[w][i][p];
                        }
                    }
                }
                T = .5;
                shader_data[66559]=T;
            }
            glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(shader_data), shader_data);
            glUniform1i(glGetUniformLocation(shaderProgram, "res"),windowres/windowScale);
            glUniform1f(glGetUniformLocation(shaderProgram, "T"),T);
            //1st pass
            float playerLocation[3];
            copypt(player.cord1,playerLocation);
            glUniform3fv(glGetUniformLocation(shaderProgram,"pl"),1,&playerLocation[0]);
            float playerLocation2[3];
            copypt(player.cord2,playerLocation2);
            glUniform3fv(glGetUniformLocation(shaderProgram,"pl2"),1,&playerLocation2[0]);
            float cameraReference[3];
            copypt(player.compass,cameraReference);
            glUniform3fv(glGetUniformLocation(shaderProgram,"camRef"),1,&cameraReference[0]);
            float duplicatePlayerLocation[3];
            copypt(player.duplicate[0].cord1,duplicatePlayerLocation);
            glUniform3fv(glGetUniformLocation(shaderProgram,"playerduploc"),1,&duplicatePlayerLocation[0]);
            float duplicatePlayerLocation2[3];
            copypt(player.duplicate[0].cord2,duplicatePlayerLocation2);
            glUniform3fv(glGetUniformLocation(shaderProgram,"playerduploc2"),1,&duplicatePlayerLocation2[0]);
            glUniform1f(glGetUniformLocation(shaderProgram, "mirrorPlayer"),player.mirror);
            glUniform1f(glGetUniformLocation(shaderProgram, "mirrorDup"),player.duplicate[0].mirror);
            glUniform1f(glGetUniformLocation(shaderProgram, "zoom"),zoom);
            glUniform1i(glGetUniformLocation(shaderProgram, "pw"),player.currentWorld);
            glUniform1i(glGetUniformLocation(shaderProgram, "duppw"),player.duplicate[0].currentWorld);
            glUniform1i(glGetUniformLocation(shaderProgram, "iopip"),player.intersectingPortalIndex);
            glUniform1i(glGetUniformLocation(shaderProgram, "diopip"),player.duplicate[0].intersectingPortalIndex);
            glBindFramebuffer(GL_FRAMEBUFFER, framebuffer1);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowres/windowScale, windowres/windowScale, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framestage1, 0);
            //std::thread render(renderFrame, windowres, windowScale, height, width, zoom, shaderProgram, framestage1, framebuffer1, postprocessfull, window, T);glfwMakeContextCurrent(window);
            //double timetest = glfwGetTime();
            glDrawArrays(GL_TRIANGLES, 0, 3);
            //timetest = glfwGetTime()-timetest;printf("%lf\n",timetest);
            //
            //2nd pass
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glUseProgram(postprocessfull);
            glUniform1i(glGetUniformLocation(postprocessfull, "height"),height);
            glUniform1i(glGetUniformLocation(postprocessfull, "width"),width);
            glUniform1i(glGetUniformLocation(postprocessfull, "windowScale"),windowScale);
            glUniform1f(glGetUniformLocation(postprocessfull, "zoom"),zoom);
            glDrawArrays(GL_TRIANGLES, 0, 3);

            // check and call events and swap the buffers
            glUseProgram(shaderProgram);//*/
        }
        glDeleteSync(fence);
        fence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
        //deltaTime=glfwGetTime()-deltaTime;
        //if(deltaTime > 1/tickRate) printf("%lf\n",1/deltaTime);

        tickCount++;
        glfwPollEvents();


        //if((tickCount&63) == 0 && glfwGetTime() > NextTickTime) printf("lag\n");
        tickCountLooped++;
        while(glfwGetTime() < NextTickTime);
        if(glfwGetTime() >= 10) {
            glfwSetTime(0);
            //printf("%lf\n",framesRendered/tickCountLooped*tickRate);
            framesRendered = 0;
            tickCountLooped = 0;
        }
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
