#include <iostream>
#include <cmath>
#define S2 1
#define E2 0
#define H2 -1

//Make changes to the input and run this
//g++ -O3 -w endpointsToCurve.cpp -o eptc && ./eptc
//the program will output a wall you can paste into the world that you can make a portal later
    //input//
    char geometry = H2;//options S2, E2, H2
    double e1[3] = {-0.0116770454, 1.91456115, 2.16001868};
    double e2[3] = {-0.994481564, 3.88759041, 4.135499};
    double r = 2;//ignored for horocycles, make negative to get the bigger of the valid arcs
    //Only expects negative inputs for arcs
    char type = 2;// 0 arc, 1 horocycle, 2 hypercycle
    /////////

double Abs(double x){
    if(x < 0) return -x;
    return x;
}

void copypt(double p1[3], double p2[3]){
    p2[0]=p1[0];p2[1]=p1[1];p2[2]=p1[2];
}

void vec3(double v[], double x, double y, double z){
    v[0]=x;v[1]=y;v[2]=z;
}

void printM(double m[][3]){
    for(int n=0;n<3;n++){
        printf("%.9lf\t%.9lf\t%.9lf\n",m[0][n],m[1][n],m[2][n]);
    }
    printf("\n");
}

void printP(double p[3]){
    printf("%.9lf, %.9lf, %.9lf\n",p[0],p[1],p[2]);
}

void matxpt(double m[][3], double p[3]){
    double A=m[0][0]*p[0]+m[1][0]*p[1]+m[2][0]*p[2],
    B=m[0][1]*p[0]+m[1][1]*p[1]+m[2][1]*p[2];
    p[2]=m[0][2]*p[0]+m[1][2]*p[1]+m[2][2]*p[2];
    p[0]=A;
    p[1]=B;
}

void mat3(double m[3][3], double a, double b, double c, double d, double e, double f, double g, double h, double i){//same format as in GLSL for ease of porting
    m[0][0]=a;m[0][1]=b;m[0][2]=c;
    m[1][0]=d;m[1][1]=e;m[1][2]=f;
    m[2][0]=g;m[2][1]=h;m[2][2]=i;
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

void backOnSphere(double p[]){
    float r=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[0]/=r;p[1]/=r;p[2]/=r;
}

double dot(double p1[3], double p2[3]){
    return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
}

void cross(double p1[3], double p2[3], double p3[3]){
    p3[0] = p1[1]*p2[2] - p1[2]*p2[1];
    p3[1] = p1[2]*p2[0] - p1[0]*p2[2];
    p3[2] = p1[0]*p2[1] - p1[1]*p2[0];
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

void lorentz(double& x, double& y, double c, double s){
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

void backOnHyperboloid(double p[]){
    p[2] = sqrt(p[0]*p[0] + p[1]*p[1] + 1);
}

void licross(double p1[3], double p2[3], double p3[3]){
    vec3(p3, p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2], p1[1]*p2[0] - p1[0]*p2[1]);
}

int main(){
    double p1[3];
    double p2[3];
    double L;
    switch(geometry){
        case 0:
        {
            e2[0] -= e1[0];e2[1] -= e1[1];
            double l = sqrt(e2[0]*e2[0] + e2[1]*e2[1]);
            p1[0] = l/2; p1[1] = -sqrt(r*r - l*l/4); p1[2] = 0;
            if(r < 0) p1[1] = -p1[1];
            copypt(p1,p2);
            p2[1] += Abs(r);
            L = (p2[0]*p2[0] + p2[1]*p2[1])/r/r;
            e2[0] /= l;e2[1] /= l;
            rotate(p1[0],p1[1],e2[0],e2[1]);
            rotate(p2[0],p2[1],e2[0],e2[1]);
            p1[0] += e1[0];p1[1] += e1[1];
            p2[0] += e1[0];p2[1] += e1[1];
        }
            break;
        case 1:
        {
            if(r < 0){
                double p[3];
                copypt(e1,p);
                copypt(e2,e1);
                copypt(p,e2);
            }
            double cosl = dot(e1,e2);
            double sinl = sqrt(1 - cosl*cosl);
            double coshalfl = sqrt((1 + cosl)/2);
            double sinhalfl = sqrt((1 - cosl)/2);
            double cosr = cos(r), sinr = sin(r);
            p1[2] = cosr/coshalfl; p1[0] = sqrt(1 - p1[2]*p1[2]); p1[1]=0;
            copypt(p1,p2);
            rotate(p2[2],p2[0],cosr,-sinr);
            L = (p2[2]*coshalfl - 1)/sinr/sinr;
            rotate(p1[2],p1[1],coshalfl,sinhalfl);
            rotate(p2[2],p2[1],coshalfl,sinhalfl);
            double rot[3][3];
            s2matto(e1,rot);
            matxpt(rot,e2);
            e2[0] /= sinl;e2[1] /= sinl;
            rotate(p1[0],p1[1],e2[1],-e2[0]);
            rotate(p2[0],p2[1],e2[1],-e2[0]);
            double roti[3][3];
            transpose(rot,roti);
            matxpt(roti,p1);
            matxpt(roti,p2);
        }
            break;
        case -1:
        {
            double coshl = lidot(e1,e2);
            double halfcoshl = sqrt((coshl + 1)/2);
            double halfsinhl = sqrt((coshl - 1)/2);
            if(type == 0) p1[2] = cosh(r);
            if(type == 1) p1[2] = 1;
            if(type == 2) p1[2] = sinh(r);
            p1[2] /= halfcoshl;
            double off =  1 - type;
            p1[0] = -sqrt(p1[2]*p1[2] - off);
            if(r < 0) p1[0] *= -1;
            p1[1] = 0;
            double k = halfcoshl*p1[2];
            if(type == 1) p2[0] = (p1[2]*p1[2] - k*k)/(2*p1[0]*k);
            else p2[0] = Abs(-k*p1[0] - p1[2]*sqrt(k*k - off));
            p2[1] = 0;
            p2[2] = sqrt(p2[0]*p2[0] + 1);
            L = p2[2]*halfcoshl;
            lorentz(p1[1],p1[2],halfcoshl,-halfsinhl);
            lorentz(p2[1],p2[2],halfcoshl,-halfsinhl);
            double iso[3][3], isoi[3][3];
            h2matto(e1,iso);
            matxpt(iso,e2);
            double sinhl = sqrt(e2[0]*e2[0] + e2[1]*e2[1]);
            e2[0] /= sinhl;e2[1] /= sinhl;
            rotate(p1[0],p1[1],-e2[1],e2[0]);
            rotate(p2[0],p2[1],-e2[1],e2[0]);
            h2invert(iso,isoi);
            matxpt(isoi,p1);
            matxpt(isoi,p2);
        }
            break;
    }
    printf("{%.9lf,%.9lf,%.9lf, %.9lf,%.9lf,%.9lf, %.9lf,%i},\n",p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], L,4+type);
}
//g++ -O3 -w endpointsToCurve.cpp -o eptc
