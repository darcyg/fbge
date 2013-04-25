/*
    Copyright 2012 thedukeofzill

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "engine.h"

//set up variables/structures
void startengine(){
    setlocale(LC_ALL,"");
    mainwin = initscr();
    keypad(mainwin, TRUE);
    nodelay(mainwin,TRUE);
    noecho();
    curs_set(0);
    start_color();
    data = make2dshortarray(1000,1000);
    zbuffer = make2ddoublearray(1000,1000);
    
    maxDepth = 10000;
    for (int j = 0; j < 1000; j++){
        for(int i = 0; i < 1000; i++){
	   zbuffer[i][j] = maxDepth;
        }
    }

    resetColors();

    stagesize = 1000; //how much we want visible in map
    dist = 100;  //eye point
    zval = 100;  //image plane
    framerate = 0;

    //initilaize color pairs
    for(int i = 0; i < 15; i++){
        for(int j = 0; j < 15; j++){
	  init_pair((i*16)+j+1,i,j);
        }
    }
}

//clear memory allocated during game
void killengine(){
    for(int i = 0; i < 1000; i++)
        free(data[i]);
    free(data);
    for(int i = 0; i < 1000; i++)
      free(zbuffer[i]);
    free(zbuffer);

    endwin();
}

//draws pixels to screen
void drawscreen(){
    for(int i = 0; i < row * 2; i += 2){
        for(int j = 0; j < col; j++){
            setcolor(data[i][j],data[i+1][j]);
            addstr("▄");
        }
    }
}

//2d vertex constructor
vertex point(double xx,double yy){
    vertex v = {xx,yy};
    return v;
}

//3d vertex constructor
vertex3 point3(double x,double y,double z){
    vertex3 v = {x,y,z};
    return v;
}

triangle tri(vertex3* v3, short color){
  triangle t;
  t.v3 = v3;
  t.color = color;
  return t;
}

poly pol(vertex* v, short cr, int s){
    poly p;
    p.color = cr;
    p.vertices = malloc(s);
    int len = p.size = s/sizeof(vertex);
    for(int i = 0; i < len; i++)
        p.vertices[i] = v[i];
    return p;
}

shape makeshape(vertex3* v3, int len, int s, short color){
    shape sh;
    sh.v2 = NULL;
    sh.v3 = v3;
    sh.size = s;
    sh.length = len;
    sh.color = color;
    return sh;
}

void swap(int* a,int* b){
    int t=*a;
    *a=*b;
    *b=t;
}

//sets color of pixel (in array, pre-draw)
void setpixelz(int rr, int cc, int color, double z){
    if(rr < r*2 && rr >= 0 && cc < c && cc >= 0 && z <= zbuffer[rr][cc]){
        data[rr][cc] = color;
        zbuffer[rr][cc] = z;
    }
}
void setpixel(int rr, int cc, int color){
    if(rr < r*2 && rr >= 0 && cc < c && cc >= 0)
        data[rr][cc] = color;
}

//draw line between two pts using midpt algorithm
void drawLine(double xx1,double yy1,double xx2,double yy2, short color){
    int dx,dy,d,incry,incre,incrne,slopegt1=0;

    int x1 = floor(xx1*c);
    int x2 = floor(xx2*c);
    int y1 = floor(yy1*r*2);
    int y2 = floor(yy2*r*2);
    
    setpixel(y1,x1,color); 

    dx=abs(x1-x2);dy=abs(y1-y2);
    if(dy>dx){
        swap(&x1,&y1);
        swap(&x2,&y2);
        swap(&dx,&dy);
        slopegt1=1;
    }
    if(x1>x2){
        swap(&x1,&x2);
        swap(&y1,&y2);
    }
    if(y1>y2)
        incry=-1;
    else
        incry=1;
    d=2*dy-dx;
    incre=2*dy;
    incrne=2*(dy-dx);
    while(x1<x2){
        if(d<=0)
            d+=incre;
        else{
            d+=incrne;
            y1+=incry;
        }
        x1++;
        if(slopegt1)
            setpixel(x1,y1,color);
        else 
            setpixel(y1,x1,color);
    }
}

//connects vertices in order, overlapping of lines may occur
void drawPoly(poly p){
    int len = p.size;
    short cr = p.color;
    vertex p1;
    vertex p2;
    for(int i = 0; i < len - 1; i++){
        p1 = p.vertices[i];
        p2 = p.vertices[i+1];
        drawLine(p1.x,p1.y,p2.x,p2.y,cr);
    }
    p1 = p.vertices[0];
    p2 = p.vertices[len-1];
    drawLine(p1.x,p1.y,p2.x,p2.y,cr);
}


short** make2dshortarray(short arraySizeX, short arraySizeY) {
    short** theArray;
    theArray = (short**) malloc(arraySizeX*sizeof(short*));
    for (int i = 0; i < arraySizeX; i++)
        theArray[i] = (short*) malloc(arraySizeY*sizeof(short));
    return theArray;
}

double** make2ddoublearray(int arraySizeX, int arraySizeY) {
  double** theArray;
  theArray = (double**) malloc(arraySizeX*sizeof(double*));
  for (int i = 0; i < arraySizeX; i++)
    theArray[i] = (double*) malloc(arraySizeY*sizeof(double));
  return theArray;
}


//resets all pixels to black
void resetColors(){
    for(int j = 0; j < col; j++){
        for(int i = 0; i < row*2; i++){
            data[i][j] = COLOR_BLACK;
        }
    }
}

//resets zbuffer
void resetZBuffer(){
  for(int j = 0; j < col; j++){
    for(int i = 0; i < row*2; i++){
      zbuffer[i][j] = maxDepth;
    }
  }
}


 
//sets colors to draw next character space (2 pixels)
void setcolor(int bg, int fg){
    attron(COLOR_PAIR((fg*16)+bg+1)); 
}

//old, draws rectangle
void drawBox(double x,double y,double l,double w,short color){
    int boxlen = floor(l*c);
    int boxwid = floor(w*r*2);
    int boxrow = floor(y*r*2-boxwid/2);
    int boxcol = floor(x*c-boxlen/2);
    if(boxrow < 0 || boxrow >= r*2) boxrow = abs(boxrow % (r*2));
    if(boxcol < 0 || boxcol >= c) boxcol = abs(boxcol % c);
    //draw horizontal
    for(int i = boxcol; i < c && i < boxcol + boxlen; i++){
        data[boxrow][i] = color;
        int temp = boxrow + boxwid;
        if(temp < r*2){
            data[temp][i] = color;
        }
    }
    ///draw vertical 
    for(int i = boxrow; i < r*2 && i < boxrow + boxwid; i++){
        data[i][boxcol] = color;
        int temp = boxcol + boxlen-1;
        if(temp < c){
            data[i][temp] = color;
        }
    } 
}

//old, draws rectangle
void drawSolidBox(double x,double y,double l,double w,short color){
    int boxlen = floor(l*c);
    int boxwid = floor(w*r*2);
    int boxrow = floor(y*r*2-boxwid/2);
    int boxcol = floor(x*c-boxlen/2);
    if(boxrow < 0 || boxrow > r*2) boxrow = abs(boxrow % (r*2));
    if(boxcol < 0 || boxcol > c) boxcol = abs(boxcol % c);
    for(int i = boxcol; i < c && i < boxcol + boxlen; i++){
        for(int j = boxrow; j <= r*2 +1&& j < boxrow + boxwid + 1; j++){
            data[j][i] = color;
        }
    } 
}


//returns new array of 3d vertices
vertex3* cubevertices(double x,double y, double z, int size){
    vertex3* v = malloc(9*sizeof(vertex3));
    int s = size/2;
    v[0] = point3(x,y,z); //center
    v[1] = point3(x-s,y-s,z-s);
    v[2] = point3(x-s,y-s,z+s);
    v[3] = point3(x+s,y-s,z+s);
    v[4] = point3(x+s,y-s,z-s);
    v[5] = point3(x-s,y+s,z-s);
    v[6] = point3(x-s,y+s,z+s);
    v[7] = point3(x+s,y+s,z+s);
    v[8] = point3(x+s,y+s,z-s);
    return v;
}

//returns cube structure
shape cube(double x,double y, double z, int size, short color){
    vertex3* cubedata = cubevertices(x,y,z,size);
    shape sh = makeshape(cubedata,9,size,color);
    return sh;
}

//determine if cube is behind camera, needs to be replaced by proper clipping
//algorithm
int visible(vertex3* array, int size){
    if(array == NULL)
        return 0;
    for(int i = 0; i < size; i++)
        if(array[i].z < -dist)
            return 0;
    return 1;
}

//find 2d coordinates of cube vertices
void translatecube(shape* sh){
    vertex3 cent = sh->v3[0];
    vertex3* newcube = 
        cubevertices(cent.x-camx,cent.y-camy,cent.z-camz,sh->size);
    rotateaxis(newcube, 9, 0.01*camrx, 0.01*camry, 0.01*camrz);
    for(int i = 0; i < 9; i++){
        //newcube[i].x += camx;
        newcube[i].y += camy;
        newcube[i].z += camz;
    }
    sh->pcenter = point3(newcube[0].x,newcube[0].y,newcube[0].z);
    if(!(visible(newcube,9))){
        vertex* v = malloc(8*sizeof(vertex));
        v[0].x = -1000000;
        sh->v2 = v;
        free(newcube);
        return;
    }
    vertex* v = malloc(8*sizeof(vertex));
    for(int i = 0; i < 8; i++){
        int denom = newcube[i+1].z + dist;
        if(!denom) denom = 1;
        v[i].x = newcube[i+1].x * (dist + zval);
        v[i].x /= denom;
        v[i].x = (v[i].x + (stagesize/2))/stagesize;

        v[i].y = newcube[i+1].y * (dist + zval);
        v[i].y /= denom;
        v[i].y = (v[i].y + (stagesize/2))/stagesize;
    }
    sh->v2 = v;
    free(newcube);
}

//draws lines in cube, wall determines whether box has Xs on sides
void drawCube(shape* sh, short wall){
  if(!wireframe){
    vertex3* cubevert = malloc(3*sizeof(vertex3));
    triangle cubetri;
        cubevert[0] = sh->v3[1];
        cubevert[1] = sh->v3[2];
        cubevert[2] = sh->v3[3];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
        cubevert[0] = sh->v3[1];
        cubevert[1] = sh->v3[3];
        cubevert[2] = sh->v3[4];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
        cubevert[0] = sh->v3[5];
        cubevert[1] = sh->v3[6];
        cubevert[2] = sh->v3[7];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);

        cubevert[0] = sh->v3[5];
        cubevert[1] = sh->v3[7];
        cubevert[2] = sh->v3[8];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
        cubevert[0] = sh->v3[1];
        cubevert[1] = sh->v3[5];
        cubevert[2] = sh->v3[8];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);


        cubevert[0] = sh->v3[1];
        cubevert[1] = sh->v3[8];
        cubevert[2] = sh->v3[4];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
        cubevert[0] = sh->v3[2];
        cubevert[1] = sh->v3[6];
        cubevert[2] = sh->v3[7];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
        cubevert[0] = sh->v3[2];
        cubevert[1] = sh->v3[7];
        cubevert[2] = sh->v3[3];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
        cubevert[0] = sh->v3[1];
        cubevert[1] = sh->v3[5];
        cubevert[2] = sh->v3[6];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);

        cubevert[0] = sh->v3[1];
        cubevert[1] = sh->v3[6];
        cubevert[2] = sh->v3[2];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
        cubevert[0] = sh->v3[4];
        cubevert[1] = sh->v3[8];
        cubevert[2] = sh->v3[7];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);

        cubevert[0] = sh->v3[4];
        cubevert[1] = sh->v3[7];
        cubevert[2] = sh->v3[3];
        cubetri = tri(cubevert,sh->color);
        drawTriangle(cubetri);
        
    free(cubevert);
    return;
  }
     

    if(sh->v2[0].x == -1000000)
        return;
    vertex* cube2d = sh->v2;
    short cl = sh->color;
    for(int i = 0; i < 4; i++){
        drawLine(cube2d[i].x,cube2d[i].y,cube2d[i+4].x,cube2d[i+4].y,cl);
        if(i < 3){
            drawLine(cube2d[i].x,cube2d[i].y,
                cube2d[i+1].x,cube2d[i+1].y,cl);
            drawLine(cube2d[i+4].x,cube2d[i+4].y,
                cube2d[i+5].x,cube2d[i+5].y,cl);
        }
        else{
            drawLine(cube2d[0].x,cube2d[0].y,cube2d[3].x,cube2d[3].y,cl);
            drawLine(cube2d[4].x,cube2d[4].y,cube2d[7].x,cube2d[7].y,cl);
        }
    }
    if(wall){
        drawLine(cube2d[0].x,cube2d[0].y,cube2d[5].x,cube2d[5].y,cl);
        drawLine(cube2d[1].x,cube2d[1].y,cube2d[6].x,cube2d[6].y,cl);
        drawLine(cube2d[2].x,cube2d[2].y,cube2d[7].x,cube2d[7].y,cl);
        drawLine(cube2d[3].x,cube2d[3].y,cube2d[4].x,cube2d[4].y,cl);
        drawLine(cube2d[0].x,cube2d[0].y,cube2d[7].x,cube2d[7].y,cl);
        drawLine(cube2d[1].x,cube2d[1].y,cube2d[4].x,cube2d[4].y,cl);
        drawLine(cube2d[2].x,cube2d[2].y,cube2d[5].x,cube2d[5].y,cl);
        drawLine(cube2d[3].x,cube2d[3].y,cube2d[6].x,cube2d[6].y,cl);
        drawLine(cube2d[0].x,cube2d[0].y,cube2d[2].x,cube2d[2].y,cl);
        drawLine(cube2d[1].x,cube2d[1].y,cube2d[3].x,cube2d[3].y,cl);
        drawLine(cube2d[4].x,cube2d[4].y,cube2d[6].x,cube2d[6].y,cl);
        drawLine(cube2d[5].x,cube2d[5].y,cube2d[7].x,cube2d[7].y,cl);
    }

}

void rotateaxis(vertex3* v, int size, double xr, double yr, double zr){
    for(int i = 0; i < size; i++){
        int x2 = v[i].x;
        int y2 = v[i].y;
        int z2 = v[i].z;

        if(zr != 0){
            v[i].x = floor(x2*cos(zr) - y2*sin(zr));
            v[i].y = floor(x2*sin(zr) + y2*cos(zr));
            if(xr != 0 || yr != 0){
                x2 = v[i].x;
                y2 = v[i].y;
            }
        }
        if(xr != 0){
            v[i].y = floor(y2*cos(xr) - z2*sin(xr));
            v[i].z = floor(y2*sin(xr) + z2*cos(xr));
            if(yr){
                y2 = v[i].y;
                z2 = v[i].z;
            }
        }
        if(yr != 0){
            v[i].z = floor(z2*cos(yr) - x2*sin(yr));
            v[i].x = floor(z2*sin(yr) + x2*cos(yr));
        }
    }
}

void rotateinplace(shape* sh, double xr, double yr, double zr){

    vertex3 cent = sh->v3[0];
    vertex3* newcube =
        cubevertices(0,0,0,sh->size);
    rotateaxis(newcube, 9, 10*xr, 10*yr, 10*zr);
    for(int i = 0; i < 9; i++){
        newcube[i].x += cent.x;
        newcube[i].y += cent.y;
        newcube[i].z += cent.z;
    }
    free(sh->v3);
    sh->v3 = newcube;
}

void drawTriangle(triangle t){
    //adjust position based on camera
    vertex3* nv = malloc(3*sizeof(vertex3));
    nv[0] = point3(t.v3[0].x-camx,t.v3[0].y-camy,t.v3[0].z-camz);
    nv[1] = point3(t.v3[1].x-camx,t.v3[1].y-camy,t.v3[1].z-camz);
    nv[2] = point3(t.v3[2].x-camx,t.v3[2].y-camy,t.v3[2].z-camz);
 
    //adjust rotation based on camera
    rotateaxis(nv, 3, 0.01*camrx, 0.01*camry, 0.01*camrz);
    for(int i = 0; i < 3; i++){
        //newcube[i].x += camx;
        nv[i].y += camy;
        nv[i].z += camz;
	if(nv[i].z < camz) return; //don't draw if behind camera
    }

    //perspective division
    vertex* v = malloc(3*sizeof(vertex));
    for(int i = 0; i < 3; i++){
        int denom = nv[i].z + dist;
        if(!denom) denom = 1;
        v[i].x = nv[i].x * (dist + zval);
        v[i].x /= denom;
        v[i].x = (v[i].x + (stagesize/2))/stagesize;

        v[i].y = nv[i].y * (dist + zval);
        v[i].y /= denom;
        v[i].y = (v[i].y + (stagesize/2))/stagesize;
    }

    //find horiz and vert bounds of triangle
    double highpt = r*2;
    double lowpt = 0;
    double rightpt = 0;
    double leftpt = c;
    double leftpty = 0;
    double topz, bottomz, leftz, rightz;
    for(int i = 0; i < 3; i++){
      double vy = v[i].y * r * 2;
      double vx = v[i].x * c;
      if (vy < highpt){
	highpt = vy;
	topz = nv[i].z;
      }
      if (vy > lowpt){
	lowpt = vy;
	bottomz = nv[i].z;
      }
      if (vx < leftpt){
        leftpt = vx;
        leftpty= vy;
	leftz = nv[i].z;
      }
      if (vx > rightpt){
        rightpt = vx;
	rightz = nv[i].z;
      }
    }

    double x0 = v[0].x * c;
    double y0 = v[0].y * r * 2;
    double x1 = v[1].x * c;
    double y1 = v[1].y * r * 2;
    double x2 = v[2].x * c;
    double y2 = v[2].y * r * 2;
    //for each y(row), find lines it intersects and points on line
    double xdepth = (rightz - leftz)/(rightpt - leftpt);
    double ydepth = (bottomz - topz)/(lowpt - highpt);
    double ztri = leftz - ydepth*(leftpty - highpt);
    for(int y = highpt; y <= lowpt; y++){
      double ztri2 = ztri;
      for(int x = leftpt; x <= rightpt; x++){
        if(tritest(x0,y0,x1,y1,x2,y2,x,y))
	  setpixelz(y,x,t.color,ztri2);
        ztri2 += xdepth;
      }
      ztri += ydepth;
    }


    //drawLine(v[0].x,v[0].y,v[1].x,v[1].y,t.color);
    //drawLine(v[0].x,v[0].y,v[2].x,v[2].y,t.color);
    //drawLine(v[2].x,v[2].y,v[1].x,v[1].y,t.color);

    free(nv);
    free(v);
  


}

int test2( double px, double py, double m, double b ) {
    if (py < m * px + b ) {
        return -1; // point is under line
    }else if ( py == m * px + b ){
        return 0; // point is on line
    } else {
        return 1; // point is over line
    }
}

int test1(double px, double py, double m,double b, double lx,double ly) {
   return (test2(px,py, m,b) == test2(lx,ly,m,b));
}

int tritest (double x0,double y0,double x1,double y1,double x2,double y2,double px, double py) {

    int line1, line2, line3;
    double m01 = (y1-y0)/(x1-x0);
    double b01 = m01 * -x1 + y1;
    double m02, m12, b02, b12;
    m02 = (y2-y0)/(x2-x0);
    m12 = (y2-y1)/(x2-x1);
    b02 = m02 * -x2 + y2;
    b12 = m12 * -x2 + y2;

    // vertical line checks

    if( x1 == x0 ) {
        line1 = ( (px <= x0) == (x2 <= x0) );
    } else {
        line1 = test1( px, py, m01, b01,x2,y2);
    }

    if( x1 == x2 ) {
        line2 = ( (px <= x2) == (x0 <= x2) );
    } else {
        line2 = test1(px,py, m12, b12,x0,y0);
    }

    if( x2 == x0 ) {
        line3 = ( (px <= x0 ) == (x1 <= x0) );} else {
        line3 = test1(px, py, m02,b02,x1,y1);
    }

    return line1 && line2 && line3;
}

