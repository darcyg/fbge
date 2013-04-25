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
#include <string.h>
#include <sys/time.h>
typedef long long int64;

static struct timeval start_time;
int64 lasttime;
int64 currtime;

void init_time() {
    gettimeofday(&start_time, NULL);
}

int64 get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (int64) (t.tv_sec - start_time.tv_sec) * 1000000
    + (t.tv_usec - start_time.tv_usec);
}



int cubecount;
int cubemax;

shape** cubes;
triangle** triangles;

//create pyramid shape
void pyramid(int x, int y, int z, int csize){
    int psize = 4;
    for(int i = x; i < x + csize*psize; i += csize){
        for(int j = y; j < y + csize*psize; j += csize){
            if(i == x || i == x + csize*(psize-1) ||
                j == y || j == y + csize*(psize-1)){
                *(cubes[cubecount]) = cube(i, z, j, csize, COLOR_RED);
                cubecount++;
            }
        }
    }

    z -= csize;
    x += csize/2;
    y += csize/2;
    psize--;
    for(int i = x; i < x + csize*psize; i += csize){
        for(int j = y; j < y + csize*psize; j += csize){
            if(i == x || i == x + csize*(psize-1) ||
                j == y || j == y + csize*(psize-1)){
                *(cubes[cubecount]) = cube(i, z, j, csize, COLOR_YELLOW);
                cubecount++;
            }
        }
    }
    z -= csize;
    x += csize/2;
    y += csize/2;
    psize--;
    for(int i = x; i < x + csize*psize; i += csize){
        for(int j = y; j < y + csize*psize; j += csize){
            if(i == x || i == x + csize*(psize-1) ||
                j == y || j == y + csize*(psize-1)){
                *(cubes[cubecount]) = cube(i, z, j, csize, COLOR_GREEN);
                cubecount++;
            }
        }
    }
    z -= csize;
    x += csize/2;
    y += csize/2;
    psize--;
    for(int i = x; i < x + csize*psize; i += csize){
        for(int j = y; j < y + csize*psize; j += csize){
            if(i == x || i == x + csize*(psize-1) ||
                j == y || j == y + csize*(psize-1)){
                *(cubes[cubecount]) = cube(i, z, j, csize, COLOR_BLUE);
                cubecount++;
            }
        }
    }
}

//create bridges
void rainbowstairs(int tracker, int cubeloc, int cubesize, int cubelocy){

        *(cubes[cubecount]) = cube(tracker, cubelocy,cubeloc+3*cubesize,
            cubesize,COLOR_RED);
        cubecount++;
        *(cubes[cubecount]) = cube(tracker, cubelocy-cubesize,cubeloc+2*cubesize,
            cubesize,COLOR_YELLOW);
        cubecount++;

        *(cubes[cubecount]) = cube(tracker, cubelocy-2*cubesize,cubeloc+cubesize,
            cubesize,COLOR_GREEN);
        cubecount++;

        *(cubes[cubecount]) = cube(tracker, cubelocy-3*cubesize,cubeloc,
            cubesize,COLOR_BLUE);
        cubecount++;

        *(cubes[cubecount]) = cube(tracker, cubelocy-2*cubesize,cubeloc-cubesize,
            cubesize,COLOR_GREEN);
        cubecount++;

        *(cubes[cubecount]) = cube(tracker, cubelocy-cubesize,cubeloc-2*cubesize,
            cubesize,COLOR_YELLOW);
        cubecount++;

        *(cubes[cubecount]) = cube(tracker, cubelocy,cubeloc-3*cubesize,
            cubesize,COLOR_RED);
        cubecount++;
}

//function to sort cubes by distance from camera
int cubesort(const void* a, const void* b){
    shape* c2 = *(shape**)a;
    shape* c1 = *(shape**)b;
    vertex3 cOne = c1->pcenter;
    vertex3 cTwo = c2->pcenter;
    if(a == NULL && b == NULL)
        return 0;
    if(b == NULL && a != NULL)
        return -1;
    if(a == NULL && b != NULL)
        return 1;
    int xa = cOne.x;
    int ya = cOne.y;
    int za = cOne.z;
    int xb = cTwo.x;
    int yb = cTwo.y;
    int zb = cTwo.z;
    //draw ground level first
    if(ya == 100 && yb != 100)
        return 1;
    if(ya != 100 && yb == 100)
        return -1;
    if(za > zb)
        return 1;
    if(za < zb)
        return -1;
    //za -= 250;
    //zb -= 250;
    double da = sqrt(xa*xa+ya*ya+za*za);
    double db = sqrt(xb*xb+yb*yb+zb*zb);
    if(da > db)
        return 1;
    if(db > da)
        return -1;
    return 0;
}


int main(int argc, char *argv[]) {
    wireframe = 1;
    cubecount = 0;
    cubemax = 10000;
    int dir = 1;
    int camdir = 1;
    short gameon = 1;
    int tracker = -100;
    int trackerprime = -tracker;
    int zdist = -50;
    int cubesize = 50;
    int cubeloc = 0;
    int cubelocy = 50;
    double stepcount = 0;
    double turncount = 0;
    camry = 0;
    camrx = 0;
    camrz = 0;
    camx = 0;
    camy = 0;
    camz = -250;
    //world objects 
    cubes = malloc(cubemax*sizeof(shape*));
    for(int i = 0; i < cubemax; i++){
        cubes[i] = malloc(sizeof(shape));
    }
    init_time();
    startengine(); //initialize stuff

    //create floor
    int ctrack = 0;
    for(int x = -1000; x < 1000; x += cubesize){
        for(int z = -1000; z < 1000; z += cubesize){
            short colr;
            if(ctrack%2 == 0) 
                colr = COLOR_BLUE; 
            else 
                colr = COLOR_GREEN;
            *(cubes[cubecount]) = 
                cube(x, 100, z, cubesize, colr);
            cubecount++;ctrack++;
        }
        ctrack++;
    }
    
    
    //place some bridges 
    for(int i = cubesize * 0; i >= 0; i -= cubesize){
        rainbowstairs(tracker - 2*i,cubeloc+cubesize*6,cubesize,cubelocy);
        rainbowstairs(-tracker + 2*i,cubeloc+cubesize*6,cubesize,cubelocy);

        rainbowstairs(tracker - 2*i,cubeloc,cubesize,cubelocy);
        rainbowstairs(-tracker + 2*i,cubeloc,cubesize,cubelocy);
        rainbowstairs(tracker - 2*i,cubeloc-cubesize*6,cubesize,cubelocy);
        rainbowstairs(-tracker + 2*i,cubeloc-cubesize*6,cubesize,cubelocy);
        rainbowstairs(tracker - 2*i,cubeloc+cubesize*12,cubesize,cubelocy);
        rainbowstairs(-tracker + 2*i,cubeloc+cubesize*12,cubesize,cubelocy);
    }
//test
    *(cubes[cubecount]) = cube(0, -50,-100, cubesize,COLOR_RED);
        cubecount++;
    *(cubes[cubecount]) = cube(0, -50,300, cubesize,COLOR_YELLOW);
        cubecount++;
    //pyramid(-500,-25,0,50);
    //variables to track movement
    int rotateleft = -1;
    int rotateright = -1;
    int movingforward = -1;
    int movingbackward = -1;
    double timepassed = 0;

    lasttime = get_time();

    
    if(argc > 1 && strcmp(argv[1],"t") == 0)
        wireframe = 0;

    //main loop
    while(gameon){
        getmaxyx(stdscr,row,col);
        r = row;
        c = col;
        
        for(int i = 0; i < cubecount; i++){
            translatecube(cubes[i]);
        }
        //sort cubes every 10 frames
        if(!(ccount%10))
            qsort((void*)cubes,cubecount,sizeof(shape*),cubesort);
        //draw cubes to buffer
        for(int i = 0; i < cubecount; i++){
            if(argc > 1 && strcmp(argv[1],"1") == 0)
                drawCube(cubes[i],1);
            else
                drawCube(cubes[i],0);
            free(cubes[i]->v2);
        }

        //copy buffer to screen
        drawscreen();

        ccount++;
        refresh();
        move(0,0);
        resetColors();
        resetZBuffer();
        currtime = get_time();
        framerate = (double)(1.0/((currtime - lasttime)*0.000001));
        timepassed = (double)(currtime - lasttime);
        lasttime = currtime;

        //get keyboard input
        ch = getch();
        if(ch != ERR){
            if(ch == KEY_RIGHT || ch =='d'){
                if(rotateright == 1 && turncount > 0){
                    turncount = 0;
                }
                else{
                    rotateleft = 1;
                    turncount = 0.15;
                }
                rotateright = -1;
            }
            else if(ch == KEY_LEFT || ch == 'a'){
                if(rotateleft == 1 && turncount > 0){
                    turncount = 0;
                }
                else{
                    rotateright = 1;
                    turncount = 0.15;
                }
                rotateleft = -1;
            }
            else if(ch == KEY_UP || ch == 'w'){
                if(movingbackward != 1){
                    movingforward = 1;
                    stepcount = 0.15;
                }
                else
                    movingbackward = -1;
            }
            else if(ch == KEY_DOWN || ch == 's'){
                if(movingforward != 1){
                    movingbackward = 1;
                    stepcount = 0.15;
                }
                else
                    movingforward = -1;
            }
            else if(ch == ' '){
                turncount = 0;
            }
            else if(ch == 'q')
                gameon = 0;
        }
        double ff = timepassed/1500;
        int speed = ff;
        if(rotateleft == 1){
            if(turncount > 0){
                camry-=ff/2;
                turncount -= 1/framerate;
            }
        }
        if(rotateright == 1){
            if(turncount > 0){
                camry+=ff/2;
                turncount -= 1/framerate;
            }

        }
        if(movingforward == 1){
            if(stepcount > 0){
                camz += speed*cos(0.01*camry);
                camx -= speed*sin(0.01*camry);
                //camz+=ff/2;
                stepcount -= 1/framerate;
            }
        }
        if(movingbackward == 1){
            if(stepcount > 0){
                camz -= speed*cos(0.01*camry);
                camx += speed*sin(0.01*camry);
                //camz-=ff/2;
                stepcount -= 1/framerate;
            }

        }
        zval = camz + 500;
        dist = -camz;
        //usleep(16000);
    }

    for(int i = 0; i < cubemax; i++){
        free(cubes[i]->v3);
        free(cubes[i]);
    }
    free(cubes);
    killengine();//prevent memory leaks
    return 0;
}
