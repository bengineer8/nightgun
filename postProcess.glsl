#version 460 core
//currently just a past through, was experimenting earlier
const float tpi=6.283185307;
const float pi=3.14159265359;
const float wallthickness=1/18.0;
uniform float zoom=1.0;
uniform int height,width;
uniform sampler2D previousstep;
out vec4 FragColor;

float brightness(vec4 color){
    return (color[0]+color[1]+color[2])*0.3;
}


void main(){
    int res=min(height,width);
    float x=gl_FragCoord.x, y=gl_FragCoord.y;
    FragColor=texture(previousstep,vec2(x/width,y/height));
    //FragColor=vec4(0,brightness(FragColor),0,1);
    /*FragColor=vec4(0,0,0,1);
    int samples=0;
    for(float X=x-1;X<x+2;X+=1){
        for(float Y=y-1;Y<y+2;Y+=1){
            if(X>=0 && X<=res && Y>=0 && Y<=res){
                FragColor+=texture(previousstep,vec2(X/width,Y/height));
                //samples+=1;
            }
        }
    }
    FragColor/=9;//*/
    //FragColor=texture(previousstep,vec2(x/width,y/height));
    //FragColor=vec4(0,brightness(FragColor),0,1);
}


/*void main(){
    int res=min(height,width);
    int pixelRad=int(res*wallthickness/(tpi*zoom)*0.5);
    float x=gl_FragCoord.x, y=gl_FragCoord.y;
    float percentNotBlack=0;
    bool nearBlack=false;
    for(float X=x-pixelRad; X<=x+pixelRad;X+=1/1){
        for(float Y=y-pixelRad;Y<=y+pixelRad;Y+=1/1){
            if(0<X && X<res && 0<Y && Y<res){
                float pixBrightness=brightness(texture(previousstep,vec2(X/width,Y/height)));
                if( (x-X)*(x-X)+(y-Y)*(y-Y)<=pixelRad*pixelRad ){
                    if(pixBrightness>0)
                        percentNotBlack+=1;
                    else
                        nearBlack=true;
                }
            }
        }
    }
    //percentNotBlack/=(zoom*zoom);
    percentNotBlack/=pi*pixelRad*pixelRad;
    percentNotBlack-=0.5;
    percentNotBlack*=1.1;
    //percentNotBlack*=percentNotBlack;percentNotBlack*=percentNotBlack;
    percentNotBlack=min(percentNotBlack,1);
    FragColor=texture(previousstep,vec2(x/width,y/height));
    if(nearBlack || x<pixelRad || y<pixelRad || x>res-pixelRad || y>res-pixelRad)
        FragColor*=vec4(percentNotBlack,percentNotBlack,percentNotBlack,1);
}*/
