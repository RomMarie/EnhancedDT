//genfeed.cpp
// Copyright © Theo Schouten and Egon van den Broek
// version 1.0, July 2013

#include <distance_transform/static/algo/libs/feed.h>

//parameters in the algorithm to obtimize execution time for a given kind of input image
// The parameters can be adjusted to obtain a minimal execution time on a given set of images.
// A way to do it is to start with an initial educated guess and then vary P1 to obtain minimal time.
// Then repeat this with P2 to P6, foolowed by going back from P5 to P1
#define SMALL1      50  // AREAbb to start filling after first bounding box, P1 in the paper
#define SMALLb     150  // AREAbb to start filling after k{1} with k<= 4, P2 in the paper
#define CUTSV      200  // cut for filling diagonal shaped A_b, P3 in the paper
#define SMALLb12v  150  // AREAbb to start filling after k{2} with k=1 search, there are v's, P4 in the paper
#define SMALLb12   200  // AREAbb to start filling after k{2} with k=1 search, there are no v's, P5 in the paper
#define SMALLnov   200  // for simple filling after the quadrant search, P6 in the paper

// Two parameters for "FILL single line" in the paper
#define CUTSY   20 // minimal size for vertical line to search k{1} with k=1 to possible reduce the line
#define CUTSX   20 // same for horizontal line

// macro's
#define AREAbb ( (maxx-minx+1)*(maxy-miny+1) ) // the size of a bounding box bb

#define MIN3(a,b,c)       ( (a) < (b) ? MIN(a,c) : MIN(b,c) )
#define MAX3(a,b,c)       ( (a) > (b) ? MAX(a,c) : MAX(b,c) )

// k(2,1) search with intersection with v's
#define UNEAR1LV \
if(v1 >= -1) {\
  m1n1=-1; if(input[y+1][x+2]==0) { m1n1=2; hasn12=1;\
   maxy=MIN(maxy,(2-2*minx)); maxx=MIN(maxx,(2-miny)/2);\
     if(v2 > -1) { j = (2*v2 + 2) / 3;  maxy = MIN(maxy, j);}\
     if(v4 > -1) { j = (v4 + 2) / 3;  maxx = MIN(maxx, j);} \
     if(v3 > -1) { j = -(2*v3 + 2); miny = MAX(miny, j); j = (v3 + 2) ;  maxx=MIN(maxx,j);} \
  } \
}else m1n1=-2;\
if(v2 >= -1) {\
  m2n1=-1; if(input[y+1][x-2]==0) { m2n1=2; hasn12=1;\
   if(m1n1 > 0){  maxy=MIN(maxy,2);} else maxy=MIN(maxy,(2+2*maxx));\
     if(v1 > -1) { j = (2*v1 + 2) / 3;  maxy = MIN(maxy, j);}\
     if(v3 > -1) { j = -(v3 + 2) / 3;  minx = MAX(minx, j);} else minx=MAX(minx,-(2-miny)/2);  \
     if(v4 > -1) { j = -(v4 + 2) ;   minx = MAX(minx, j); j = -(2*v4 + 2) ; miny = MAX(miny, j);   }  }  \
} else m2n1=-2;\
if(v3 >= -1) {\
  m3n1=-1; if(input[y-1][x-2]==0) { m3n1=2; hasn12=1;\
   if(m2n1 > 0){  minx=MAX(minx,-1);} else minx=MAX(minx,-(2+maxy)/2);\
      if(v4 > -1) { j = -(2*v4 + 2) / 3;  miny = MAX(miny, j);} else miny=MAX(miny,-(2+2*maxx));\
      if(v2 > -1) { j = -(v2 + 2) / 3;  minx = MAX(minx, j);} \
      if(v1 > -1) { j = (2*v1 + 2) ;   maxy=MIN(maxy,j); j = -(v1 + 2) ;  minx = MAX(minx, j);} \
  } \
} else m3n1=-2;\
if(v4 >= -1) {\
  m4n1=-1; if(input[y-1][x+2]==0) { m4n1=2; hasn12=1;\
   if(m3n1 > 0){  miny=MAX(miny,-2);} else miny=MAX(miny,-(2-2*minx));\
   if(m1n1 > 0){  maxx=MIN(maxx,1);} else  maxx=MIN(maxx,(2+maxy)/2);\
     if(v3 > -1) { j = -(2*v3 + 2) / 3;  miny = MAX(miny, j);}\
     if(v1 > -1) { j =  (v1 + 2) / 3;  maxx = MIN(maxx, j);} \
     if(v2 > -1) { j = (v2 + 2) ;  maxx=MIN(maxx,j); j = (2*v2 + 2) ; maxy=MIN(maxy,j); } \
  }\
} else m4n1=-2;

// k(2,1) search no intersection with v's
#define UNEAR1 \
if(v1 >= -1) {\
  m1n1=-1; if(input[y+1][x+2]==0) { m1n1=2; hasn12=1;\
   maxy=MIN(maxy,(2-2*minx)); maxx=MIN(maxx,(2-miny)/2);\
  } \
}else m1n1=-2;\
if(v2 >= -1) {\
  m2n1=-1; if(input[y+1][x-2]==0) { m2n1=2; hasn12=1;\
    minx=MAX(minx,-(2-miny)/2); \
   if(m1n1 > 0){  maxy=MIN(maxy,2);} else maxy=MIN(maxy,(2+2*maxx));\
  }\
} else m2n1=-2;\
if(v3 >= -1) {\
  m3n1=-1; if(input[y-1][x-2]==0) { m3n1=2; hasn12=1;\
   miny=MAX(miny,-(2+2*maxx)); \
   if(m2n1 > 0){  minx=MAX(minx,-1);} else minx=MAX(minx,-(2+maxy)/2);\
  } \
} else m3n1=-2;\
if(v4 >= -1) {\
  m4n1=-1; if(input[y-1][x+2]==0) { m4n1=2; hasn12=1;\
   \
   if(m3n1 > 0){  miny=MAX(miny,-2);} else miny=MAX(miny,-(2-2*minx));\
   if(m1n1 > 0){  maxx=MIN(maxx,1);} else  maxx=MIN(maxx,(2+maxy)/2);\
  }\
} else m4n1=-2;

// k(1,2) search with intersection with v's
#define UNEAR2LV \
if(v1 >= -1) {\
  m1n2=-1; if(input[y+2][x+1]==0) { m1n2=2; hasn12=1;\
   maxy=MIN(maxy,(2-minx)/2); maxx=MIN(maxx,2-2*miny);\
     if(v2 > -1) { j = (v2 + 2) / 3;  maxy = MIN(maxy, j);}\
     if(v4 > -1) { j = (2*v4 + 2) / 3;  maxx = MIN(maxx, j);} \
     if(v3 > -1) { j = (v3 + 2) ;   maxy=MIN(maxy,j); j = -(2*v3 + 2) ;  minx = MAX(minx, j);}\
  }  \
}else m1n2=-2;\
if(v2 >= -1) {\
  m2n2=-1; if(input[y+2][x-1]==0) { m2n2=2;  hasn12=1;\
  minx=MAX(minx,-(2-2*miny));\
   if(m1n2 > 0){  maxy=MIN(maxy,1);} else maxy=MIN(maxy,(2+maxx)/2);\
     if(v1 > -1) { j = (v1 + 2) / 3;  maxy = MIN(maxy, j);}\
     if(v3 > -1) { j = -(2*v3 + 2) / 3;  minx = MAX(minx, j);} \
     if(v4 > -1) { j = (v4 + 2) ;   maxy=MIN(maxy,j);  j = (2*v4 + 2)  ;   maxx=MIN(maxx,j);}\
  } \
} else m2n2=-2;\
if(v3 >= -1) {\
  m3n2=-1; if(input[y-2][x-1]==0) { m3n2=2; hasn12=1; \
   miny=MAX(miny,-(2+maxx)/2); \
   if(m2n2 > 0){  minx=MAX(minx,-2);} else minx=MAX(minx,-(2+2*maxy));\
      if(v4 > -1) { j = -(v4 + 2) / 3;  miny = MAX(miny, j);}\
      if(v2 > -1) { j = -(2*v2 + 2) / 3;  minx = MAX(minx, j);} \
      if(v1 > -1) { j = -(v1 + 2) ;   miny = MAX(miny, j);  j = 2*v1+2;   maxx=MIN(maxx,j);}\
  }\
} else m3n2=-2;\
if(v4 >= -1) {\
  m4n2=-1; if(input[y-2][x+1]==0) { m4n2=2;  hasn12=1;\
   if(m3n2 > 0){  miny=MAX(miny,-1);} else miny=MAX(miny,-(2-minx)/2);\
   if(m1n2 > 0){  maxx=MIN(maxx,2);} else maxx=MIN(maxx,2+2*maxy);\
     if(v3 > -1) { j = -(v3 + 2) / 3;  miny = MAX(miny, j);}\
     if(v1 > -1) { j =  (2*v1 + 2) / 3;  maxx = MIN(maxx, j);} \
     if(v2 > -1) { j = -(v2 + 2) ;   miny = MAX(miny, j);  j = -(2*v2 + 2) ;   minx = MAX(minx, j); }\
  }\
} else m4n2=-2;

// k(1,2) search, no intersection with v's
#define UNEAR2 \
if(v1 >= -1) {\
  m1n2=-1; if(input[y+2][x+1]==0) { m1n2=2; hasn12=1;\
   maxy=MIN(maxy,(2-minx)/2); maxx=MIN(maxx,2-2*miny);\
  }  \
}else m1n2=-2;\
if(v2 >= -1) {\
  m2n2=-1; if(input[y+2][x-1]==0) { m2n2=2;  hasn12=1;\
  minx=MAX(minx,-(2-2*miny));\
   if(m1n2 > 0){  maxy=MIN(maxy,1);} else maxy=MIN(maxy,(2+maxx)/2);\
  } \
} else m2n2=-2;\
if(v3 >= -1) {\
  m3n2=-1; if(input[y-2][x-1]==0) { m3n2=2; hasn12=1; \
   miny=MAX(miny,-(2+maxx)/2); \
   if(m2n2 > 0){  minx=MAX(minx,-2);} else minx=MAX(minx,-(2+2*maxy));\
  }\
} else m3n2=-2;\
if(v4 >= -1) {\
  m4n2=-1; if(input[y-2][x+1]==0) { m4n2=2;  hasn12=1;\
   if(m3n2 > 0){  miny=MAX(miny,-1);} else miny=MAX(miny,-(2-minx)/2);\
   if(m1n2 > 0){  maxx=MIN(maxx,2);} else maxx=MIN(maxx,2+2*maxy);\
  }\
} else m4n2=-2;



int Cfeed::genfeed(const CImg<uchar>& inpcimg, CImg<int>& output)
{
 const  int xsize = inpcimg.width, ysize = inpcimg.height; // row (x) and column (y) size of image
 const  int xsizem1 = xsize-1, ysizem1 = ysize -1;

 int x, y; // used for looping over the image
 int lx, ly; // usually used for locally scanning the image around an (x,y) point
 int j; // temporary local variable
 int low, high; // used for filling the area around a border pixel, low and high point on a row
 int ed, edy; // used in filling

 int minx, maxx, miny, maxy;  // the extend of a bounding box
 int v1, v2, v3, v4; int hasv; // for the 45 degree k(1,1) search lines in the 4 quadrants
 // vq = -2 : a bisection line can not cut into the bb
// vq = -1 : a bisection line was not yet found
// vq >= 0 : the value of the right side of eq[7] in the paper, defining the bisection line x+y=v1 etc
 int m1n1, m2n1, m3n1, m4n1; int m1n2, m2n2, m3n2, m4n2; int hasn12; // mqn1(2) for the k(2,1) resp k(1,2) search lines in the four quadrants
// when mqni > 0 it defines a bisection line

 int xnext, xprev; // next and previous border pixel on a row
 int h2; //used is search for the top of a bb

// define 2 matrices, input[y][x] giving access to the input image and dt[y][x] to the output distance image
 uchar* input[ysize];
    for(y=0;y<ysize;y++)input[y] = (uchar*)inpcimg.data +y*xsize;
 int * dt[ysize];
    for(y=0;y<ysize;y++)dt[y] = ( int *)output.data +y*xsize;

 { // initialize the output distances: 0 for an object pixel, otherwise a value larger than the maximal distance
    const  int MAXVAL = ysize*ysize + xsize*xsize;
    uchar *in= ( uchar *)inpcimg.data;  int *dtp = ( int *)output.data; //&dt[0][0];
    for(lx=0; lx < ysize*xsize; lx++) if(*in++ == 0)*dtp++=0; else *dtp++ = MAXVAL;
 }

 int  ylast[xsize], ynext[xsize]; // auxiliary vectors for the quadrant search
 for(x=0; x < xsize; x++) ylast[x] = 0; // the previous border pixel in a column
 for(x=0; x < xsize; x++) ynext[x] = ysize; // possible the next border pixel in a column


for ( y = 4; y < ysize-4; y++ ) { // loop over the rows

 for(x=4; x < xsize-4 ; x++) if( input[y][x]==0) break; // find first object pixel

 xprev= -x;

 for ( ; x < xsize-4; x=xnext ) { // always at a border
  minx = xprev;
  if(input[y-1][x] == 0)  {miny=0;} // we reached a top border pixel
  else {
      if(ylast[x]==0)
          miny=-y;
      else miny= - (y-ylast[x]-1)/2 ;
  } // determine bottum of bb from ylast[x]
  ylast[x] =y;

  if( input[y+1][x] ==0) {h2=y+1;  maxy=0;}  // short search for top of bb
  else if(input[y+2][x]==0) {h2=y+2;  maxy=1;}
  else if(input[y+3][x]==0) {h2=y+3;  maxy=1;}
  else if(input[y+4][x]==0) {h2=y+4;  maxy=2;}
  else { h2 = ysize; maxy = ysizem1-y; } // not found, set to maximum
  ynext[x]=h2;

// find right side of the bb; also find next border pixel on the row
  if( input[y][x+1]==0) {
        maxx=0;
        lx=x; do {lx++;} while ( input[y][lx+1]==0 && input[y-1][lx] ==0 && input[y+1][lx] ==0); // move onto next border pixel
        xnext=lx; xprev=0;
  }else {
      maxx=xsizem1-x;
      for( ly=1 ; ly < xsize-4-x; ly++) if(input[y][x+ly] == 0 ) { maxx=ly/2;break;}
      xnext = x+ly; xprev= - (ly-1)/2;
  }

// bb is single vertical line
  if(minx==0&&maxx==0) { // faster to test on v's when part of line is big enough
    if(maxy >  CUTSY && (input[y+1][x-1]==0 || input[y+1][x+1]==0) ) dt[y+1][x]=1;
    else for( edy=1,ly=1 ; ly<=maxy; ly++) { if( edy < dt[y+ly][x] ) dt[y+ly][x] = edy; else break; edy+=2*ly+1; }
    if(miny < -CUTSY && (input[y-1][x-1]==0 || input[y-1][x+1]==0) ) dt[y-1][x]=1;
    else for(edy=1,ly=-1; ly>=miny; ly--) { if( edy < dt[y+ly][x] ) dt[y+ly][x] = edy; else break; edy+=1-2*ly;}
    continue;
  }
//bb is single horizontal line
  if(miny==0&&maxy==0) {
    if(maxx >  CUTSX && (input[y+1][x+1]==0 || input[y-1][x+1]==0) ) dt[y][x+1]=1;
    else for(edy=1,ly=1 ; ly<=maxx; ly++) { if( edy < dt[y][x+ly] ) dt[y][x+ly] = edy; else break; edy+=2*ly+1; }
    if(minx < -CUTSX && (input[y+1][x-1]==0 || input[y-1][x-1]==0) ) dt[y][x-1]=1;
    else for(edy=1,ly=-1; ly>=minx; ly--) { if( edy < dt[y][x+ly] ) dt[y][x+ly] = edy; else break; edy+=1-2*ly;}
    continue;
  }


  if(h2==ysize && AREAbb > SMALL1) { // continue search for top of bb
      uchar *p = &input[y+5][x];
      for(  ly=5 ; ly < ysizem1-y-4; ly++, p += xsize) if(*p == 0 ) { maxy=ly/2; ynext[x]=y+ly; break;}
  }

  if( AREAbb < SMALL1) { // P1 in the paper
      edy = (minx-1)*(minx-1);
      int *p = &dt[y+miny][x];
      for(ly=miny; ly<= maxy; ly++) {
          ed= edy+ly*ly;
          for(lx=minx; lx <=maxx; lx++) { ed += 2*lx-1; if(ed <p[lx]) p[lx] = ed;}
          p += xsize;
      }
      continue;
  }

// search k{1} with k<=4
  hasv=0;
  if( maxy+maxx >1) {
    if(     input[y+1][x+1] == 0 ) v1=1; else if(input[y+2][x+2] == 0 ) v1=2;
    else if(input[y+3][x+3] == 0 ) v1=3; else if(input[y+4][x+4] == 0 ) v1=4; else v1= -1;
    if(v1 > 0) {hasv=1; maxy=MIN(maxy,v1-minx);maxx=MIN(maxx,v1-miny);}
  } else v1=-2;
  if(maxy-minx > 1) {
    if(     input[y+1][x-1] == 0 ) v2=1; else if(input[y+2][x-2] == 0 ) v2=2;
    else if(input[y+3][x-3] == 0 ) v2=3; else if(input[y+4][x-4] == 0 ) v2=4; else v2= -1;
    if(v2 > 0) {hasv=1; maxy=MIN(maxy,v2+maxx);minx=MAX(minx,-(v2-miny));
        if(v1 > 0){ j=(v1+v2)/2; maxy=MIN(maxy,j); } }
  } else v2=-2;
  if(-miny-minx > 0) {
    if(     input[y-1][x-1] == 0 ) v3=0; else if(input[y-2][x-2] == 0 ) v3=1;
    else if(input[y-3][x-3] == 0 ) v3=2; else if(input[y-4][x-4] == 0 ) v3=3; else v3= -1;
    if(v3 >= 0) {hasv=1; miny=MAX(miny,-(v3+maxx));minx=MAX(minx, -(v3+maxy));
        if(v2 > 0){ j=(v2+v3)/2; minx=MAX(minx,-j);} }
  } else v3=-2;
  if(maxx-miny > 0) {
    if(     input[y-1][x+1] == 0 ) v4=0; else if(input[y-2][x+2] == 0 ) v4=1;
    else if(input[y-3][x+3] == 0 ) v4=2; else if(input[y-4][x+4] == 0 ) v4=3; else v4= -1;
    if(v4 >= 0) {hasv=1; miny=MAX(miny,-(v4-minx));maxx=MIN(maxx,v4+maxy);
        if(v3 >= 0){ j=(v3+v4)/2; miny=MAX(miny,-j);}
        if(v1 > 0) { j=(v1+v4)/2; maxx=MIN(maxx,j);} }
  } else v4=-2;


  if( hasv>0  ) { // one or more bisection lines from k{1}found

    if(AREAbb < SMALLb ) { // small enough bb
        if(minx==0&&maxx==0){
            for(ly=miny; ly<=maxy; ly++) {edy=ly*ly; if( edy < dt[ly+y][x] ) dt[ly+y][x] = edy; }
        } else  {
            edy = (minx-1)*(minx-1);
            int *p = &dt[y+miny][x];
            for(ly=miny; ly<= maxy; ly++) {
                ed= edy+ly*ly;
                for(lx=minx; lx <=maxx; lx++) { ed += 2*lx-1; if(ed <p[lx]) p[lx] = ed;}
                p += xsize;
            }
        }
        continue;
    }


    if(v2==1 && v4==0  && (maxx-minx) < CUTSV) { // small enough diagonal A_b
        low=MAX(minx,miny); high=MIN(maxx,maxy);
        ed=2*low*low; edy=ed+1 -2*low;
        if(low-1 >= minx) {if(edy <dt[y+low][x+low-1]) dt[y+low][x+low-1] =edy;}
        if(ed <dt[y+low][x+low]) dt[y+low][x+low] =ed;
        for(ly=low+1; ly<= high; ly++) {
            ed+= 4*ly-2;  if(ed <dt[y+ly][x+ly]) dt[y+ly][x+ly] =ed;
            edy+= 4*ly-4; if(edy <dt[y+ly][x+ly-1]) dt[y+ly][x+ly-1] =edy; // slower: edy=ed-2*ly+1;
        }
        if(high+1 <= maxy){ edy=high*high + (high+1)*(high+1); if(edy <dt[y+high+1][x+high]) dt[y+high+1][x+high] =edy;}
        continue;
    }

    if(v1==1 && v3==0 && (maxx-minx) < CUTSV) {
        low=MAX(-maxx,miny); high=MIN(-minx,maxy);
        if((-low+1) <= maxx) {edy= low*low+(low-1)*(low-1); if(edy <dt[y+low][x+(-low+1)]) dt[y+low][x+(-low+1)] =edy; }
        ed=2*low*low; if(ed <dt[y+low][x-low]) dt[y+low][x-low] =ed;
        ed=2*low*low; edy=ed+1 -2*low;
        for(ly=low+1; ly<= high; ly++) {
            ed+= 4*ly-2; if(ed <dt[y+ly][x-ly]) dt[y+ly][x-ly] =ed;
            edy+= 4*ly-4; if(edy <dt[y+ly][x+(-ly+1)]) dt[y+ly][x+(-ly+1)] =edy;
        }
        if(high+1 <= maxy){ edy=high*high +(high+1)*(high+1); if(edy <dt[y+high+1][x-high]) dt[y+high+1][x-high] =edy;}
        continue;
    }

    hasn12=0; UNEAR1LV; UNEAR2LV; // search k{2} with k=1, include intersection with v's

    if(hasn12==1 && AREAbb < SMALLb12v ) { // small enough bb
        if(minx==0&&maxx==0){
            for(ly=miny; ly<=maxy; ly++) {edy=ly*ly; if( edy < dt[ly+y][x] ) dt[ly+y][x] = edy; }
        } else  {
            int *p = &dt[y+miny][x]; edy= (minx-1)*(minx-1);
            for(ly=miny; ly<= maxy; ly++) {
                ed= edy+ly*ly;
                for(lx=minx; lx <=maxx; lx++) { ed += 2*lx-1; if(ed <p[lx]) p[lx] = ed;}
                p += xsize;
            }
        }
        continue;
    }


  } // end of test on any v's
  else { // there are no v's

    hasn12=0;  UNEAR1; UNEAR2; // seach lines k{2} with k=1, no intersection test with v (they are not there)

    if(hasn12==1 && AREAbb < SMALLb12 ) { // small enough bb
        if(minx==0&&maxx==0){
            for(ly=miny; ly<=maxy; ly++) {edy=ly*ly; if( edy < dt[ly+y][x] ) dt[ly+y][x] = edy; }
        } else  {
            edy = (minx-1)*(minx-1);  int *p = &dt[y+miny][x];
            for(ly=miny; ly<= maxy; ly++) {
                ed= edy+ly*ly;
                for(lx=minx; lx <=maxx; lx++) { ed += 2*lx-1; if(ed < p[lx] ) p[lx] = ed;}
                p += xsize;
            }
        }
        continue;
    }


  } // end of split between with or without v's

//start quadrant search

// an object point at (lx,ly) from the borer point gives a bisection line equation:
//  ly*y + lx*x = (lx^2+ly^2 + ud) /2
//   with ud=0 in quadrants 1 and 2, and -1 in quadrants 3 and 4
// instead of keeping per quadrant the lowest right hand size, we keep an approxiation |lx|+|ly|
// variables used per quadrant: cai=ly, cb1i=lx, cabsi=|lx|+|ly|

#define CHIGH 1000000000 // value of cabsi indicating no bisection line found

  int ca1=0,cb1=0,ca2=0,cb2=0,ca3=0,cb3=0,ca4=0,cb4=0, cabs1=CHIGH, cabs2=CHIGH, cabs3=CHIGH, cabs4=CHIGH;

  if(v1==-1) { // only do the quadrant search if no v has been found
    if(m1n1==2 ) { ca1= 1; cb1= 2; cabs1=3;} // take a k{2} bisection line
    else if(m1n2==2 ) { ca1= 2; cb1= 1; cabs1=3;}
    else  { // the real search
        if(ynext[x+1]< ysize) {
            ly = ynext[x+1]-y;
            if(ly > 0) { cabs1=1+ly;  ca1=ly; cb1=1; }
        }
        const int ix1= MIN(xsize-5-x, MAX(2*maxx,maxx+maxy) );
        int ed=MIN(ix1,cabs1); // indicates the last lx to be searched
        for(lx=2; lx < ed ; lx++) if( ynext[x+lx]< ysize) {
            ly = ynext[x+lx]-y;
            if(ly>0 && lx+ly < cabs1) {cabs1=lx+ly; ed = MIN(ix1,cabs1); ca1=ly; cb1=lx; }
        }
    }
  }

  if( v2==-1) {
    if(m2n1==2) {ca2= 1; cb2=-2; cabs2=3;}
    else if(m2n2==2) {ca2= 2; cb2=-1; cabs2=3;}
    else  {
        if(ynext[x-1]< ysize ) { ly = ynext[x-1]-y; if(ly>0) { cabs2=1+ly;  ca2=ly; cb2=-1; } }
        const int ix2= MIN( x-5, -MIN(2*minx,minx-maxy) ); int ed= MIN(ix2,cabs2);
        for(lx=2; lx < ed; lx++) if( ynext[x-lx]< ysize ) {
            ly = ynext[x-lx]-y; if(ly>0 && lx+ly < cabs2) {cabs2=lx+ly; ed= MIN(ix2,cabs2); ca2=ly; cb2=-lx; }
        }
    }
  }

  if(v3==-1) {
    if(m3n1==2) {ca3=-1; cb3=-2; cabs3=3;}
    else if(m3n2==2) {ca3=-2; cb3=-1; cabs3=3;}
    else  {
        if(ylast[x-1] > 0 ) { ly = ylast[x-1]-y; if(ly<0) { cabs3=1-ly;  ca3=ly; cb3=-1; } }
        const int ix3=MIN(x-5, -MIN(2*minx,minx+miny) ); int ed=MIN(ix3,cabs3);
        for(lx=2; lx < ed; lx++) if( ylast[x-lx] > 0 ) {
            ly = ylast[x-lx]-y; if(ly<0 && lx-ly < cabs3) {cabs3=lx-ly; ed=MIN(ix3,cabs3); ca3=ly; cb3=-lx; }
        }
    }
  }

  if(v4==-1){
    if(m4n1==2) {ca4=-1; cb4= 2; cabs4=3;}
    else if(m4n2==2) {ca4=-2; cb4= 1; cabs4=3;}
    else  {
        if (ylast[x+1] > 0) { ly = ylast[x+1]-y; { cabs4=1-ly;  ca4=ly; cb4=1; } }
        const int ix4= MIN(xsize-5-x, MAX(2*maxx,maxx-miny) ); int ed=MIN(ix4,cabs4);
        for(lx=2; lx < ed; lx++) if( ylast[x+lx] > 0) {
            ly = ylast[x+lx]-y; if(lx-ly < cabs4) {cabs4=lx-ly;ed=MIN(ix4,cabs4);  ca4=ly; cb4=lx; }
        }
    }
  }

// define cri as (lx^2+ly^2 + ud) /2 in there is a bisection line in quadrant i, otherwise CHIGH
const int cr1= (cabs1 < CHIGH) ? (ca1*ca1+cb1*cb1)/2 :CHIGH;
const int cr2= (cabs2 < CHIGH) ? (ca2*ca2+cb2*cb2)/2 :CHIGH;
const int cr3= (cabs3 < CHIGH) ? (ca3*ca3+cb3*cb3-1)/2 :CHIGH;
const int cr4= (cabs4 < CHIGH) ? (ca4*ca4+cb4*cb4-1)/2 :CHIGH;

// reducing bb
  if(cabs1 < CHIGH) {
    if(v2>-1) { if ( (v2*cb1+cr1) < maxy*(cb1+ca1) ) maxy= (v2*cb1+cr1) / (cb1+ca1); }
    else if(cabs2 < CHIGH) { if ( (cr2*cb1-cr1*cb2) < maxy*(ca2*cb1-ca1*cb2) ) maxy= (cr2*cb1-cr1*cb2) / (ca2*cb1-ca1*cb2); }
    if(v4>-1) { if ( (v4*ca1+cr1) < maxx*(cb1+ca1) ) maxx= (v4*ca1+cr1) / (cb1+ca1); }
    else if(cabs4 < CHIGH) { if ( (cr4*ca1-cr1*ca4) < maxx*(ca1*cb4-ca4*cb1) ) maxx= (cr4*ca1-cr1*ca4) / (ca1*cb4-ca4*cb1); }
  } else if( v1 > -1) {
    if(cabs2 < CHIGH) { if ( (cr2-v1*cb2) < maxy*(ca2-cb2) ) maxy= (cr2-v1*cb2) / (ca2-cb2); }
    if(cabs4 < CHIGH) { if ( (cr4-v1*ca4) < maxx*(cb4-ca4) ) maxx= (cr4-v1*ca4) / (cb4-ca4); }
  }

  if(cabs3 < CHIGH) {
    if(v2>-1) { if ( (cr2+v2*ca2) < minx*(ca3+cb3) ) minx= (cr2+v2*ca2) / (ca3+cb3); }
    else if(cabs2 < CHIGH) { if ( (cr2*ca3-cr3*ca2) > minx*(ca3*cb2-ca2*cb3) ) minx= (cr2*ca3-cr3*ca2) / (ca3*cb2-ca2*cb3); }
    if(v4>-1) { if ( (cr3-v4*cb3) < miny*(ca3+cb3) ) miny= (cr3-v4*cb3) / (ca3+cb3); }
    else if(cabs4 < CHIGH) { if ( (cr4*cb3-cr3*cb4) > miny*(ca4*cb3-ca3*cb4) ) miny= (cr4*cb3-cr3*cb4) / (ca4*cb3-ca3*cb4); }
  } else if (v3>-1) {
    if(cabs2 < CHIGH) { if ( (cr2+v3*ca2) < minx*(cb2-ca2) ) minx= (cr2+v3*ca2) / (cb2-ca2); }
    if(cabs4 < CHIGH) { if ( (cr4+v3*cb4) < miny*(ca4-cb4) ) miny= (cr4+v3*cb4) / (ca4-cb4); }
  }


  if( AREAbb < SMALLnov ) { // bb small enough, fill it completely
    edy = (minx-1)*(minx-1);
    int *p = &dt[y+miny][x];
    for(ly=miny; ly<= maxy; ly++) {
        ed= edy+ly*ly;
        for(lx=minx; lx <=maxx; lx++) { ed += 2*lx-1; if(ed < p[lx] ) p[lx] = ed;}
        p += xsize;
    }
  }
  else  { // fill A_b
     int *p = &dt[y+miny][x];
     edy = miny*miny;
     for(ly=miny; ly<= maxy; ly++) {
        if(v1 > -1) high = MIN(maxx,v1-ly); else if(cabs1 < CHIGH) high= MIN(maxx, (cr1-ca1*ly)/cb1); else high=maxx;
        if(v2 > -1) low = MAX(minx,ly-v2);  else if(cabs2 < CHIGH) low = MAX(minx, (cr2-ca2*ly)/cb2); else low=minx;
        if(v3 > -1) low = MAX(low,-ly-v3); else if(cabs3 < CHIGH) low = MAX(low, (cr3-ca3*ly)/cb3);
        if(v4 > -1) high = MIN(high,v4+ly); else if(cabs4 < CHIGH) high= MIN(high,(cr4-ca4*ly)/cb4);
        ed = edy +(low-1)*(low-1); edy += 2*ly+1;
        for(lx=low; lx <=high; lx++) { ed += 2*lx-1; if(ed < p[lx]) p[lx] = ed;}
        p += xsize;
    }
  }

 } // end loop over x
} // end of loop over y,


return 0;
}
