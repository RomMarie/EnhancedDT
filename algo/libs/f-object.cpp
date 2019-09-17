// fobject.cpp
// Copyright © Theo Schouten and Egon van den Broek
// version 1.0, July 2013

#include <distance_transform/static/algo/libs/feed.h>

// the papameter of the algoritm

#define CUTB   36    // cut on initial bb, P1 in the paper
#define CUTSV 125    // cut on A_b after small v, P2 in the paper
#define CUTS  125    // cut bb after small v, P3 in the paper
#define CUT12 2500   // P5 in the paper
#define CUT12L 10000 // P7 in the paper

#define D11 12  // P4 in the paper
#define D12 8   // P6 in the paper
#define D14 4   // P8 in the paper
#define D34 4   // also P8 in the paper, might be set different to D14

// usefull macro's
#define AREAbb ( (maxx-minx+1)*(maxy-miny+1) ) // the area on the bounding box

#define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#define MIN(a,b) (( (a) < (b) ) ? (a) : (b))
#define MIN3(a,b,c)       ( (a) < (b) ? MIN(a,c) : MIN(b,c) )


//for searching k{1} with k=5 step D11 (parameter P5 in paper)
#define DISTV12  if(v1 == -1){ ly=MIN3(maxx+maxy,(xsizem1-x), (ysizem1-y)); \
    for( lx=5; lx <= ly; lx += D11)if(input[y+lx][x+lx]==0) { while(input[y+lx-1][x+lx-1]==0)lx--; v1=lx;   \
     maxy=MIN(maxy,v1-minx);  maxx=MIN(maxx,v1-miny);\
     if(m2n2>0) { j = (v1 + m2n2) / 3;  maxy = MIN(maxy, j);}\
     else if(v2 >= 0){ j=(v1+v2)/2; maxy=MIN(maxy,j);}\
     else if(m2n1>0) { j = (2*v1 + m2n1) / 3;  maxy = MIN(maxy, j);}\
     if(m4n1>0) { j =  (v1 + m4n1) / 3;  maxx = MIN(maxx, j);} \
     else if(v4 >= 0){ j=(v1+v4)/2; maxx=MIN(maxx,j);}\
     else if(m4n2>0){ j =  (2*v1 + m4n2) / 3;  maxx = MIN(maxx, j);} \
     break;   } }\
   if(v2 == -1){ ly=MIN3(maxy-minx, x, (ysizem1-y));\
    for( lx=5; lx <= ly; lx += D11)if(input[y+lx][x-lx]==0) {while(input[y+lx-1][x-lx+1]==0)lx--; v2=lx;   \
     maxy=MIN(maxy,v2+maxx); minx=MAX(minx,-(v2-miny));\
     if(m1n2>0){ j = (v2 + m1n2) / 3;  maxy = MIN(maxy, j);}\
     else if(v1 >= 0){ j=(v1+v2)/2; maxy=MIN(maxy,j);}\
     else if(m1n1>0) { j = (2*v2 + m1n1) / 3;  maxy = MIN(maxy, j);}\
     if(m3n1>0) { j = -(v2 + m3n1) / 3;  minx = MAX(minx, j);} \
     else if(v3 >= 0){ j=(v2+v3)/2; minx=MAX(minx,-j);}\
     else if(m3n2>0){ j = -(2*v2 + m3n2) / 3;  minx = MAX(minx, j);} \
     break;   } }\
   if(v3 == -1){ ly=MIN3(-minx-miny,x, y);\
    for( lx=5; lx <= ly; lx += D11)if(input[y-lx][x-lx]==0) {while(input[y-lx+1][x-lx+1]==0)lx--; v3=lx-1;  \
     miny=MAX(miny,-(v3+maxx)); minx=MAX(minx,-(v3+maxy));\
     if(m2n1>0) { j = -(v3 + m2n1) / 3;  minx = MAX(minx, j);} \
     else if(v2 >= 0){ j=(v2+v3)/2; minx=MAX(minx,-j);}\
     else if(m2n2>0) { j = -(2*v3 + m2n2) / 3;  minx = MAX(minx, j);} \
     if(m4n2>0) { j = -(v3 + m4n2) / 3;  miny = MAX(miny, j);}\
     else if(v4 >= 0){ j=(v3+v4)/2; miny=MAX(miny,-j);}\
     else if(m4n1>0) { j = -(2*v3 + m4n1) / 3;  miny = MAX(miny, j);}\
     break;   } }\
   if(v4 == -1){ly=MIN3(maxx-miny,(xsizem1-x), y);\
    for( lx=5; lx <= ly; lx += D11)if(input[y-lx][x+lx]==0) {while(input[y-lx+1][x+lx-1]==0)lx--; v4=lx-1;  \
     miny=MAX(miny,-(v4-minx)); maxx=MIN(maxx,v4+maxy);\
     if(m3n2>0) { j = -(v4 + m3n2) / 3;  miny = MAX(miny, j);}\
     else if(v3 >= 0){ j=(v3+v4)/2; miny=MAX(miny,-j);}\
     else if(m3n1>0) { j = -(2*v4 + m3n1) / 3;  miny = MAX(miny, j);}\
     if(m1n1>0) { j = (v4 + m1n1) / 3;  maxx = MIN(maxx, j);} \
     else if(v1 >= 0){ j=(v1+v4)/2; maxx=MIN(maxx,j);}\
     else if(m1n2>0) { j = (2*v4 + m1n2) / 3;  maxx = MIN(maxx, j);} \
     break;   } }

// for searching k{4} with k=1
#define UNEARV(i) \
if(maxx+maxy > 1) {  if(input[y+sy##i][x+sx##i]==0) { m1n##i=mc##i/2;\
     maxy=MIN(maxy,(m1n##i-mb##i*minx)/ma##i); maxx=MIN(maxx,(m1n##i-ma##i*miny)/mb##i);\
     if(v2 > -1) { j = (mb##i*v2 + m1n##i) / (ma##i+mb##i);  maxy = MIN(maxy, j);}\
     if(v4 > -1) { j = (ma##i*v4 + m1n##i) / (ma##i+mb##i);  maxx = MIN(maxx, j);} \
     if(v3 > -1) { j = (mb##i*v3 + m1n##i) / (ma##i-mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                   j = (ma##i*v3 + m1n##i) / (-ma##i+mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
  }else m1n##i=-1;} else m1n##i=-2;\
if(maxy-miny > 1){  if(input[y+sy##i][x-sx##i]==0) { m2n##i=mc##i/2;\
     maxy=MIN(maxy,(m2n##i+mb##i*maxx)/ma##i); minx=MAX(minx,-((m2n##i-ma##i*miny)/mb##i));\
     if(m1n##i >= 0){ j=(m1n##i+m2n##i)/(2*ma##i); maxy=MIN(maxy,j);} \
     if(v1 > -1) { j = (mb##i*v1 + m2n##i) / (ma##i+mb##i);  maxy = MIN(maxy, j);}\
     if(v3 > -1) { j = -(ma##i*v3 + m2n##i) / (ma##i+mb##i);  minx = MAX(minx, j);} \
     if(v4 > -1) { j = (mb##i*v4 + m2n##i) / (ma##i-mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                   j = (ma##i*v4 + m2n##i) / (ma##i-mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
  }else m2n##i=-1;} else m2n##i=-2;\
if(miny+minx < 1){  if(input[y-sy##i][x-sx##i]==0) { m3n##i=(mc##i-1)/2;\
     miny=MAX(miny,-((m3n##i+mb##i*maxx)/ma##i)); minx=MAX(minx,-((m3n##i+ma##i*maxy)/mb##i));\
     if(m2n##i >= 0){ j=(m2n##i+m3n##i)/(2*mb##i); minx=MAX(minx,-j);} \
      if(v4 > -1) { j = -(mb##i*v4 + m3n##i) / (ma##i+mb##i);  miny = MAX(miny, j);}\
      if(v2 > -1) { j = -(ma##i*v2 + m3n##i) / (ma##i+mb##i);  minx = MAX(minx, j);} \
      if(v1 > -1) { j = (mb##i*v1 + m3n##i) / (-ma##i+mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                    j = (ma##i*v1 + m3n##i) / (ma##i-mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
  }else m3n##i=-1;} else m3n##i=-2;\
if(maxx-miny > 1){  if(input[y-sy##i][x+sx##i]==0) { m4n##i=(mc##i-1)/2;\
     miny=MAX(miny,-((m4n##i-mb##i*minx)/ma##i)); maxx=MIN(maxx,(m4n##i+ma##i*maxy)/mb##i);\
     if(m3n##i >= 0){ j=(m3n##i+m4n##i)/(2*ma##i); miny=MAX(miny,-j);}\
     if(m1n##i >= 0){ j=(m1n##i+m4n##i)/(2*mb##i); maxx=MIN(maxx,j);}\
     if(v3 > -1) { j = -(mb##i*v3 + m4n##i) / (ma##i+mb##i);  miny = MAX(miny, j);}\
     if(v1 > -1) { j =  (ma##i*v1 + m4n##i) / (ma##i+mb##i);  maxx = MIN(maxx, j);} \
     if(v2 > -1) { j = (mb##i*v2 + m4n##i) / (-ma##i+mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                   j = (ma##i*v2 + m4n##i) / (-ma##i+mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
  }else m4n##i=-1;} else m4n##i=-2;

// for searching k(2,1) with k=1,2
#define UNEAR12LV1 \
if(v1 >= -1) {\
  m1n1=-1;; if(input[y+1][x+2]==0) m1n1=2; else if(input[y+2][x+4]==0) m1n1=5;\
  if(m1n1>0) { \
   maxy=MIN(maxy,(m1n1-2*minx)); maxx=MIN(maxx,(m1n1-miny)/2);\
     if(v2 > -1) { j = (2*v2 + m1n1) / 3;  maxy = MIN(maxy, j);}\
     if(v4 > -1) { j = (v4 + m1n1) / 3;  maxx = MIN(maxx, j);} \
     if(v3 > -1) { j = -(2*v3 + m1n1); miny = MAX(miny, j); j = (v3 + m1n1) ;  maxx=MIN(maxx,j);} \
  } \
}else m1n1=-2;\
if(v2 >= -1) {\
  m2n1=-1; if(input[y+1][x-2]==0) m2n1=2; else if(input[y+2][x-4]==0) m2n1=5;\
  if(m2n1>0){\
   maxy=MIN(maxy,(m2n1+2*maxx)); minx=MAX(minx,-((m2n1-miny)/2));\
   if(m1n1 >= 0){ j=(m1n1+m2n1)/(2); maxy=MIN(maxy,j);}\
     if(v1 > -1) { j = (2*v1 + m2n1) / 3;  maxy = MIN(maxy, j);}\
     if(v3 > -1) { j = -(v3 + m2n1) / 3;  minx = MAX(minx, j);} \
     if(v4 > -1) { j = -(v4 + m2n1) ;   minx = MAX(minx, j); j = -(2*v4 + m2n1) ; miny = MAX(miny, j);   }  }  \
} else m2n1=-2;\
if(v3 >= -1) {\
  m3n1=-1; if(input[y-1][x-2]==0) m3n1=2; else if(input[y-2][x-4]==0) m3n1=4;\
  if(m3n1>0){ \
   miny=MAX(miny,-((m3n1+2*maxx))); minx=MAX(minx,-((m3n1+maxy)/2));\
   if(m2n1 >= 0){ j=(m2n1+m3n1)/(4); minx=MAX(minx,-j);}\
      if(v4 > -1) { j = -(2*v4 + m3n1) / 3;  miny = MAX(miny, j);}\
      if(v2 > -1) { j = -(v2 + m3n1) / 3;  minx = MAX(minx, j);} \
      if(v1 > -1) { j = (2*v1 + m3n1) ;   maxy=MIN(maxy,j); j = -(v1 + m3n1) ;  minx = MAX(minx, j);} \
  } \
} else m3n1=-2;\
if(v4 >= -1) {\
  m4n1=-1; if(input[y-1][x+2]==0) m4n1=2; else if(input[y-2][x+4]==0) m4n1=4;\
  if(m4n1>0){\
   miny=MAX(miny,-((m4n1-2*minx))); maxx=MIN(maxx,(m4n1+maxy)/2);\
   if(m3n1 >= 0){ j=(m3n1+m4n1)/(2); miny=MAX(miny,-j);}\
   if(m1n1 >= 0){ j=(m1n1+m4n1)/(4); maxx=MIN(maxx,j);}\
     if(v3 > -1) { j = -(2*v3 + m4n1) / 3;  miny = MAX(miny, j);}\
     if(v1 > -1) { j =  (v1 + m4n1) / 3;  maxx = MIN(maxx, j);} \
     if(v2 > -1) { j = (v2 + m4n1) ;  maxx=MIN(maxx,j); j = (2*v2 + m4n1) ; maxy=MIN(maxy,j); } \
  }\
} else m4n1=-2;

// for searching k(1,2) with k=1,2
#define UNEAR12LV2 \
if(v1 >= -1) {\
  m1n2=-1; if(input[y+2][x+1]==0) m1n2=2; else if(input[y+4][x+2]==0) m1n2=5;\
  if(m1n2>0) {\
   maxy=MIN(maxy,(m1n2-minx)/2); maxx=MIN(maxx,(m1n2-2*miny));\
     if(v2 > -1) { j = (v2 + m1n2) / 3;  maxy = MIN(maxy, j);}\
     if(v4 > -1) { j = (2*v4 + m1n2) / 3;  maxx = MIN(maxx, j);} \
     if(v3 > -1) { j = (v3 + m1n2) ;   maxy=MIN(maxy,j); j = -(2*v3 + m1n2) ;  minx = MAX(minx, j);}\
  }  \
}else m1n2=-2;\
if(v2 >= -1) {\
  m2n2=-1; if(input[y+2][x-1]==0) m2n2=2; else if(input[y+4][x-2]==0) m2n2=5;\
  if(m2n2>0){ \
   maxy=MIN(maxy,(m2n2+maxx)/2); minx=MAX(minx,-((m2n2-2*miny)));\
   if(m1n2 >= 0){ j=(m1n2+m2n2)/(4); maxy=MIN(maxy,j);}\
     if(v1 > -1) { j = (v1 + m2n2) / 3;  maxy = MIN(maxy, j);}\
     if(v3 > -1) { j = -(2*v3 + m2n2) / 3;  minx = MAX(minx, j);} \
     if(v4 > -1) { j = (v4 + m2n2) ;   maxy=MIN(maxy,j);  j = (2*v4 + m2n2)  ;   maxx=MIN(maxx,j);}\
  } \
} else m2n2=-2;\
if(v3 >= -1) {\
  m3n2=-1; if(input[y-2][x-1]==0) m3n2=2; else if(input[y-4][x-2]==0) m3n2=4;\
  if(m3n2>0){ \
   miny=MAX(miny,-((m3n2+maxx)/2)); minx=MAX(minx,-((m3n2+2*maxy)));\
   if(m2n2 >= 0){ j=(m2n2+m3n2)/(2); minx=MAX(minx,-j);}\
      if(v4 > -1) { j = -(v4 + m3n2) / 3;  miny = MAX(miny, j);}\
      if(v2 > -1) { j = -(2*v2 + m3n2) / 3;  minx = MAX(minx, j);} \
      if(v1 > -1) { j = -(v1 + m3n2) ;   miny = MAX(miny, j);  j = 2*v1+m3n2;   maxx=MIN(maxx,j);}\
  }\
} else m3n2=-2;\
if(v4 >= -1) {\
  m4n2=-1; if(input[y-2][x+1]==0) m4n2=2; else if(input[y-4][x+2]==0) m4n2=4;\
  if(m4n2>0){ \
   miny=MAX(miny,-((m4n2-minx)/2)); maxx=MIN(maxx,(m4n2+2*maxy));\
   if(m3n2 >= 0){ j=(m3n2+m4n2)/(4); miny=MAX(miny,-j);}\
   if(m1n2 >= 0){ j=(m1n2+m4n2)/(2); maxx=MIN(maxx,j);}\
     if(v3 > -1) { j = -(v3 + m4n2) / 3;  miny = MAX(miny, j);}\
     if(v1 > -1) { j =  (2*v1 + m4n2) / 3;  maxx = MIN(maxx, j);} \
     if(v2 > -1) { j = -(v2 + m4n2) ;   miny = MAX(miny, j);  j = -(2*v2 + m4n2) ;   minx = MAX(minx, j); }\
  }\
} else m4n2=-2;

//for searching k{2} with k=3 step D12 (parameter P6 in paper) and k{4} with k=2 step D14 and D34 (parameter P8 in paper)
#define DISTOV(i)  if(m1n##i == -1){ j=x+MIN3( (ma##i*maxy + mb##i*maxx)*2*sx##i/mc##i,(xsizem1-x), (sx##i*(ysizem1-y))/sy##i);\
    for( ly=y+sty##i, lx=x+stx##i; lx <= j; ly += sdy##i, lx += sdx##i)if(input[ly][lx]==0) { while(input[ly-sy##i][lx-sx##i]==0)lx-=sx##i,ly-=sy##i;  m1n##i=(mc##i*(lx-x)/sx##i)/2;  \
     maxy=MIN(maxy,(m1n##i-mb##i*minx)/ma##i); maxx=MIN(maxx,(m1n##i-ma##i*miny)/mb##i);\
     if(v2 > -1) { j = (mb##i*v2 + m1n##i) / (ma##i+mb##i);  maxy = MIN(maxy, j);}\
     if(v4 > -1) { j = (ma##i*v4 + m1n##i) / (ma##i+mb##i);  maxx = MIN(maxx, j);} \
     if(v3 > -1) { j = (mb##i*v3 + m1n##i) / (ma##i-mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                   j = (ma##i*v3 + m1n##i) / (-ma##i+mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
     if(m2n##i >= 0){ j=(m1n##i+m2n##i)/(2*ma##i); maxy=MIN(maxy,j);}\
     if(m4n##i >= 0){ j=(m1n##i+m4n##i)/(2*mb##i); maxx=MIN(maxx,j);}\
     break;   } }\
   if(m2n##i == -1){ j = x-MIN3((ma##i*maxy - mb##i*minx)*2*sx##i/mc##i, x, (sx##i*(ysizem1-y))/sy##i);\
    for( ly=y+sty##i, lx=x-stx##i; lx >= j; ly += sdy##i, lx -= sdx##i)if(input[ly][lx]==0) { while(input[ly-sy##i][lx+sx##i]==0)lx+=sx##i,ly-=sy##i;m2n##i=(mc##i*(x-lx)/sx##i)/2;  \
     maxy=MIN(maxy,(m2n##i+mb##i*maxx)/ma##i); minx=MAX(minx,-((m2n##i-ma##i*miny)/mb##i));\
     if(v1 > -1) { j = (mb##i*v1 + m2n##i) / (ma##i+mb##i);  maxy = MIN(maxy, j);}\
     if(v3 > -1) { j = -(ma##i*v3 + m2n##i) / (ma##i+mb##i);  minx = MAX(minx, j);} \
     if(v4 > -1) { j = (mb##i*v4 + m2n##i) / (ma##i-mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                   j = (ma##i*v4 + m2n##i) / (ma##i-mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
     if(m1n##i >= 0){ j=(m1n##i+m2n##i)/(2*ma##i); maxy=MIN(maxy,j);}\
     if(m3n##i >= 0){ j=(m2n##i+m3n##i)/(2*mb##i); minx=MAX(minx,-j);}\
     break;   } }\
   if(m3n##i == -1){ j = MIN3( (-ma##i*miny - mb##i*minx)*2*sx##i/mc##i,x, (sx##i*y)/sy##i);\
    for( ly=sty##i, lx=stx##i; lx <= j; ly += sdy##i, lx += sdx##i)if(input[y-ly][x-lx]==0) { while(input[y-ly+sy##i][x-lx+sx##i]==0)lx-=sx##i,ly-=sy##i;m3n##i=(mc##i*(lx)/sx##i-1)/2;\
     miny=MAX(miny,-((m3n##i+mb##i*maxx)/ma##i)); minx=MAX(minx,-((m3n##i+ma##i*maxy)/mb##i));\
      if(v4 > -1) { j = -(mb##i*v4 + m3n##i) / (ma##i+mb##i);  miny = MAX(miny, j);}\
      if(v2 > -1) { j = -(ma##i*v2 + m3n##i) / (ma##i+mb##i);  minx = MAX(minx, j);} \
      if(v1 > -1) { j = (mb##i*v1 + m3n##i) / (-ma##i+mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                    j = (ma##i*v1 + m3n##i) / (ma##i-mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
     if(m2n##i >= 0){ j=(m2n##i+m3n##i)/(2*mb##i); minx=MAX(minx,-j);}\
     if(m4n##i >= 0){ j=(m3n##i+m4n##i)/(2*ma##i); miny=MAX(miny,-j);}\
     break;   } }\
   if(m4n##i == -1){ j = MIN3( (-ma##i*miny + mb##i*maxx)*2*sx##i/mc##i,(xsizem1-x), (sx##i*y)/sy##i);\
     for( ly=sty##i, lx=stx##i; lx <= j; ly += sdy##i, lx += sdx##i)if(input[y-ly][x+lx]==0) { while(input[y-ly+sy##i][x+lx-sx##i]==0)lx-=sx##i,ly-=sy##i;m4n##i=(mc##i*lx/sx##i-1)/2;\
     miny=MAX(miny,-((m4n##i-mb##i*minx)/ma##i)); maxx=MIN(maxx,(m4n##i+ma##i*maxy)/mb##i);\
     if(v3 > -1) { j = -(mb##i*v3 + m4n##i) / (ma##i+mb##i);  miny = MAX(miny, j);}\
     if(v1 > -1) { j =  (ma##i*v1 + m4n##i) / (ma##i+mb##i);  maxx = MIN(maxx, j);} \
     if(v2 > -1) { j = (mb##i*v2 + m4n##i) / (-ma##i+mb##i);  if (j < 0) miny = MAX(miny, j); else maxy=MIN(maxy,j); \
                   j = (ma##i*v2 + m4n##i) / (-ma##i+mb##i);  if (j < 0) minx = MAX(minx, j); else maxx=MIN(maxx,j);} \
     if(m3n##i >= 0){ j=(m3n##i+m4n##i)/(2*ma##i); miny=MAX(miny,-j);}\
     if(m1n##i >= 0){ j=(m1n##i+m4n##i)/(2*mb##i); maxx=MIN(maxx,j);}\
     break;   } }

// using vq's to define a range of pixels to be fed on a row
#define VCUTF \
	if(v1 > -1) high = MIN(maxx,v1-ly); else high=maxx;\
	if(v2 > -1) low = MAX(minx,ly-v2); else low=minx;\
	if(v3 > -1) low = MAX(low,-(ly+v3));\
	if(v4 > -1) high = MIN(high,v4+ly);

// using mqni's to define a range of pixels to be fed on a row
#define OCUT(i) \
	if(m1n##i > -1) high = MIN(high,(m1n##i-ma##i*ly)/mb##i);\
	if(m2n##i > -1) low = MAX(low,(ma##i*ly-m2n##i)/mb##i);\
	if(m3n##i > -1) low = MAX(low,(-ma##i*ly-m3n##i)/mb##i);\
	if(m4n##i > -1) high = MIN(high,(m4n##i+ma##i*ly)/mb##i);

// check whether found vq's and mqni's still cut into the bb
#define VCHECK \
if(v1 >= (maxx+maxy))v1=-2; if(v2>=(maxy-minx))v2=-2; if(v3>= (-minx-miny))v3=-2; if(v4 >= (maxx-miny))v4=-2;
#define OCHECK(i) \
	if(m1n##i >=  mb##i*maxx+ma##i*maxy) m1n##i= -2;	if(m2n##i >= -mb##i*minx+ma##i*maxy) m2n##i= -2;\
	if(m3n##i >= -mb##i*minx-ma##i*miny) m3n##i= -2;	if(m4n##i >=  mb##i*maxx-ma##i*miny) m4n##i= -2;

// fill diagonal shaped (2 lines) A_b when v2=1 and v4=0
#define FILLv2v4 \
	low=MAX(minx,miny); high=MIN(maxx,maxy);\
	ed=0; edy=1;\
	for(ly=1; ly<= high; ly++) {\
	    edy+= 4*ly-4; if(edy <dt[y+ly][x+ly-1]) dt[y+ly][x+ly-1] =edy;\
	    ed+= 4*ly-2;  if(ed <dt[y+ly][x+ly]) dt[y+ly][x+ly] =ed; else {\
	        for(ly++;ly<= high; ly++) {edy+= 4*ly-4; if(edy <dt[y+ly][x+ly-1]) dt[y+ly][x+ly-1] =edy; else high=MIN(high,ly);  } }\
	}\
	if(high+1 <= maxy){ edy=high*high + (high+1)*(high+1); if(edy <dt[y+high+1][x+high]) dt[y+high+1][x+high] =edy;}\
	ed=0; edy=1;\
	for(ly=-1; ly > low; ly--) {\
	    edy+= -4*ly; if(edy <dt[y+ly][x+ly-1]) dt[y+ly][x+ly-1] =edy;\
	    ed+= -4*ly-2;  if(ed <dt[y+ly][x+ly]) dt[y+ly][x+ly] =ed; else {\
	        for(ly--; ly > low; ly--) { edy+= -4*ly; if(edy <dt[y+ly][x+ly-1]) dt[y+ly][x+ly-1] =edy; else low=MAX(low,ly);} }\
	}\
	if(low-1 >= minx) {edy=low*low+(low-1)*(low-1); if(edy <dt[y+low][x+low-1]) dt[y+low][x+low-1] =edy;}\
	           ed=2*low*low; if(ed <dt[y+low][x+low]) dt[y+low][x+low] =ed;

// fill diagonal shaped (2 lines) A_b when v1=1 and v3=0
#define FILLv1v3\
	low=MAX(-maxx,miny); high=MIN(-minx,maxy);\
    ed=0; edy=1;\
	for(ly=1; ly <= high; ly++) {\
		edy+= 4*ly-4; if(edy <dt[y+ly][x+(-ly+1)]) dt[y+ly][x+(-ly+1)] =edy;\
		ed+= 4*ly-2; if(ed <dt[y+ly][x-ly]) dt[y+ly][x-ly] =ed;\
		  else { for(ly++; ly <= high; ly++) { edy+= 4*ly-4; if(edy <dt[y+ly][x+(-ly+1)]) dt[y+ly][x+(-ly+1)] =edy; else high=MIN(high,ly); } }\
	}\
	if(high+1 <= maxy){ edy=high*high +(high+1)*(high+1); if(edy <dt[y+high+1][x-high]) dt[y+high+1][x-high] =edy;}\
	ed=0; edy=1;\
	for(ly=-1; ly > low; ly--) {\
		edy+= -4*ly; if(edy <dt[y+ly][x+(-ly+1)]) dt[y+ly][x+(-ly+1)] =edy;\
		ed+= -4*ly-2; if(ed <dt[y+ly][x-ly]) dt[y+ly][x-ly] =ed;\
		  else { for(ly--; ly > low; ly--) { edy+= -4*ly; if(edy <dt[y+ly][x+(-ly+1)]) dt[y+ly][x+(-ly+1)] =edy; else low=MAX(low,ly); }  }\
	}\
	if((-low+1) <= maxx) {edy= low*low+(low-1)*(low-1); if(edy <dt[y+low][x+(-low+1)]) dt[y+low][x+(-low+1)] =edy; }\
	                ed=2*low*low; if(ed <dt[y+low][x-low]) dt[y+low][x-low] =ed;


int Cfeed::fobject(const CImg<uchar>& inpcimg, CImg<int>& output)
{
    int xsize = inpcimg.width, ysize = inpcimg.height; // x and y of image
    int xsizem1 = xsize-1, ysizem1 = ysize -1;

const int MAXVAL = ysize*ysize + xsize*xsize; // larger than the maximal squared ed in the image

int x,y; // usually the (x,y) of the current border point
int bi; // index for looping over the borders of a row
int lx,ly; // coordinates of a point relative to the current border pixel
int j; // local variable
int low, high; // usually the low and high point on a row to fill
int minx, maxx, miny, maxy; // bounding box (bb) size
int ed, edy; // used for calculatting the square ed

// search lines constants and variables
// an object point at (k*sx,k*sy) from the considered border point
// gives a bisection line:
// ma*y +mb*x = (mc * k +ud)/2  with ma = sy, mb = sx, mc=ma^2+mb^2
// ud=0 in quadrants 1 and 2, and -1 in quadrants 3 and 4

int v1, v2, v3, v4; // for the 45 degree k(1,1) search lines in the 4 quadrants
// vq = -2 : a bisection line can not cut into the bb
// vq = -1 : a bisection line was not yet found
// vq >= 0 : the value of (mc * k +ud)/2
// note that ma, mb are 1 and mc=2

// search lines k{2}
const int ma1=1, mb1=2, mc1=5, sx1=2, sy1=1;
const int ma2=2, mb2=1, mc2=5, sx2=1, sy2=2;
int m1n1, m2n1, m3n1, m4n1; // mqn1 as vq for k(2,1)
int m1n2, m2n2, m3n2, m4n2; // mqn2 as vq for k(1,2)
const int stx1=6, sty1=3, sdx1=2*D12, sdy1=D12 ; // start point and steps for the long distance search for k(2,1)
const int stx2=3, sty2=6, sdx2=D12, sdy2=2*D12 ;

// search lines k{4}
const int ma3=1, mb3=4, mc3=17, sx3=4, sy3=1; //k(4,1)
const int ma4=4, mb4=1, mc4=17, sx4=1, sy4=4; //k(1,4)
int m1n3, m2n3, m3n3, m4n3;
int m1n4, m2n4, m3n4, m4n4;
const int stx3=8, sty3=2, sdx3=4*D14, sdy3=D14 ;
const int stx4=2, sty4=8, sdx4=D14, sdy4=4*D14 ;

const int ma5=3, mb5=4, mc5=25, sx5=4, sy5=3; //k(4,3)
const int ma6=4, mb6=3, mc6=25, sx6=3, sy6=4; //k(3,4)
int m1n5, m2n5, m3n5, m4n5;
int m1n6, m2n6, m3n6, m4n6;
const int stx5=8, sty5=6, sdx5=4*D34, sdy5=3*D34 ;
const int stx6=6, sty6=8, sdx6=3*D34, sdy6=4*D34 ;



// construct input matrix pointing to inpcimg input image
uchar * input[ysize];
for(y=0;y<ysize;y++)input[y] = (uchar*)inpcimg.data +y*xsize;

//construct dt matrix, receiving the minimum distance, use directly the space of the output image:
int* dt[ysize];
for(y=0;y<ysize;y++)dt[y] = ( int *)output.data +y*xsize;

// start constructing the list of borders and initializing the output matrix
struct bor {int x; int minx; int maxx; int miny; int maxy;}; // border structure
int bn[ysize];   //number of borders per y row

{ // put code in separate block
    int *dtp= &dt[0][0]; // pointer for initializing the output matrix with the ssquared distance to the closest object pixel in the row
    bor* bdp= (bor *)smdata; // pointer to list of border stored in smdata, a buffer defined in feed.h
    int yprev[xsize];  for(x=0; x < xsize; x++) yprev[x] = 0; // previous border pixel in a column
    bor *nfill[xsize]; // pointer to border where maxy still must be filled in
    int xprev;
    int up, down; // pixel above and below the current border pixel

    std::cout<<"REUSSI"<<std::endl;
    y=0;while( y< 4*xsize ) { *dtp++= MAXVAL; y++;} // initialize first 4 rows with no object pixels
    for(y=4;y<ysize-4;y++) { //loop over the rows
        for(x=4; x < xsize-4; x++)if(input[y][x]==0) break; // find first object pixel on the row
        if(x >= xsize-4) {bn[y]=0; x=0;while(x<xsize) {*dtp++ = MAXVAL;x++;} continue;} // no object pixel on the row
        for(lx=x, ed=x*x; lx>=0;lx--){*dtp++ = ed; ed+=1-2*lx;}// fill output upto this first border pixel
        j=0;  bdp->minx=-x;
        xprev=x; up=input[y+1][xprev]; down=input[y-1][xprev];
        x++; while(x<xsize-4) { // continue loop over a row
         if(input[y][x]==0) { // now at the next border pixel in the row
            maxx=(x-xprev)/2; // bb maxx for the previous border pixel

            if( down+ up>0) {
                if(down > 0) {
                    if(yprev[xprev]==0) {bdp->miny = -y;}
                    else { nfill[xprev]->maxy= ( y-yprev[xprev] )/2; bdp->miny = -( y-yprev[xprev] -1)/2;}
                } else { bdp->miny=1;}
                if(up > 0) { yprev[xprev]=y; nfill[xprev]=bdp;} else { bdp->maxy=-1;}
                bdp->maxx= maxx;bdp->x=xprev; j++;bdp++;
            } // do not put the border in the list if up and down are both object pixels

            if(maxx==0) { // previous pixel is the previous border pixel
                while(input[y][x+1]==0 && input[y-1][x]==0 && input[y+1][x]==0) {x++,*dtp++ = 0;} // skip untill next border, set x to that
                *dtp++= 0;
                bdp->minx=0;
             }else { // fill output from previous border pixel to current one
                  for(lx=ed=1;lx<=maxx;lx++) {*dtp++ = ed; ed+=2*lx+1;}
                  bdp->minx =lx=-( x-xprev-1)/2; for(ed=lx*lx;lx<=0;lx++) {*dtp++ = ed; ed+=2*lx+1;} // fills the border pixel
             }
             xprev=x; up=input[y+1][xprev]; down=input[y-1][xprev];
          }
          x++;
        } // continue loop over the row
        for(lx=ed=1;lx<xsize-xprev;lx++) {*dtp++ = ed; ed+=2*lx+1;} // initialize output from last border pixel to end of the row
        if( down+ up>0) { // handle last border pixel on the row
            if(down > 0) {
                if(yprev[xprev]==0) {bdp->miny = -y;}
                else { nfill[xprev]->maxy= ( y-yprev[xprev] )/2; bdp->miny = -( y-yprev[xprev] -1)/2;}
            } else { bdp->miny=1;}
            if(up > 0) { yprev[xprev]=y; nfill[xprev]=bdp;} else { bdp->maxy=-1;}
            bdp->maxx= (xsizem1-xprev);bdp->x=xprev; j++;bdp++;
        }
        bn[y]=j; //number of border pixels on row y
    } // end of loop over the rows
    x=0; while(x < xsize) { if(yprev[x]>0) nfill[x]->maxy= ysizem1-yprev[x]; x++;} // fill in maxy for last border pixel in a columnn
    y=0; while(y<4*xsize){*dtp++= MAXVAL;y++;} // initialize last 4 rows with no object pixels
}

std::cout<<"REUSSI"<<std::endl;
// start looping over the borders
bor *bdp= (bor *)smdata; // pointer to the list of border pixels
for ( y = 4; y < ysize-4; y++ ) if(bn[y]>0){ //loop over rows with border pixels
 for(bi=0;bi<bn[y];bi++) {
  x= bdp->x;  minx = bdp->minx; maxx= bdp->maxx; miny = bdp->miny; maxy= bdp++->maxy;

//check for small enough bb
  if(AREAbb < CUTB) { // bb is small enough, fill it
      edy = (minx-1)*(minx-1);
      for(ly=miny; ly<= maxy; ly++) {
          ed= edy+ly*ly;
          int *p = &dt[ly+y][minx+x];
          for(lx=minx; lx <=maxx; lx++,p++) { ed += 2*lx-1; if(ed <*p) *p = ed;}
      }
      continue;
  }

  if(minx==0&&maxx==0) { // bb is part of single vertical line, fill it
      if(maxy>0 && dt[y+1][x] > 1) // test on v1=1 or v2=1
        for(ly=edy=1 ; ly<=maxy; ly++) { if( edy < dt[y+ly][x] ) dt[y+ly][x] = edy; else break; edy+=2*ly+1; }
      //if(miny<0 && dt[y-1][x]>1) // test on v3=0 or v4=0; is slower because the loop below breaks out early
        for(ly=edy=1; ly<= -miny; ly++) { if( edy < dt[y-ly][x] ) dt[y-ly][x] = edy; else break; edy+=2*ly+1;}
      continue;
  }

  if(miny>=0 && maxy <= 0) continue; // bb is part of single horizontal line, is allready filled

  if( x==4 || x==xsize-5 || y==4 || y==ysize-5) { // no need to search further because the edge of the image does not contain object pixels
    if(minx==0 && maxx==0){
        for(ly=miny; ly<=maxy; ly++) {edy=ly*ly; if( edy < dt[ly+y][x] ) dt[ly+y][x] = edy;}
    } else  {
        edy = minx*minx;
        for(ly=miny; ly<= maxy; ly++) {
            ed= edy+ly*ly;
            int *p = &dt[ly+y][minx+x];
            for(lx=minx; lx <=maxx; lx++,p++) { if(ed <*p) *p = ed; ed += 2*lx+1; }
        }
    }
    continue;
  }

//search k{1} with k<=4

  if( maxy+maxx >1) {
    if(     input[y+1][x+1] == 0 ) v1=1; else if(input[y+2][x+2] == 0 ) v1=2;
    else if(input[y+3][x+3] == 0 ) v1=3; else if(input[y+4][x+4] == 0 ) v1=4; else v1= -1;
    if(v1 > -1) {
        maxy=MIN(maxy,v1-minx);maxx=MIN(maxx,v1-miny);} //change bb using v1
    }  else v1=-2;
  if(maxy-minx > 1) {
    if(     input[y+1][x-1] == 0 ) v2=1; else if(input[y+2][x-2] == 0 ) v2=2;
    else if(input[y+3][x-3] == 0 ) v2=3; else if(input[y+4][x-4] == 0 ) v2=4; else v2= -1;
    if(v2 > -1) {
        maxy=MIN(maxy,v2+maxx);minx=MAX(minx,-(v2-miny));
        if(v1 >= 0){ j=(v1+v2)/2; maxy=MIN(maxy,j); } } // use intersection of v1 and v2 to check for reducing the bb
  } else v2=-2;
  if(-miny-minx > 0) {
    if(     input[y-1][x-1] == 0 ) v3=0; else if(input[y-2][x-2] == 0 ) v3=1;
    else if(input[y-3][x-3] == 0 ) v3=2; else if(input[y-4][x-4] == 0 ) v3=3; else v3= -1;
    if(v3 > -1) {
        miny=MAX(miny,-(v3+maxx));minx=MAX(minx, -(v3+maxy));
        if(v2 >= 0){ j=(v2+v3)/2; minx=MAX(minx,-j);} }
  } else v3=-2;
  if(maxx-miny > 0) {
    if(     input[y-1][x+1] == 0 ) v4=0; else if(input[y-2][x+2] == 0 ) v4=1;
    else if(input[y-3][x+3] == 0 ) v4=2; else if(input[y-4][x+4] == 0 ) v4=3; else v4= -1;
    if(v4 > -1) {
        miny=MAX(miny,-(v4-minx));maxx=MIN(maxx,v4+maxy);
        if(v3 >= 0){ j=(v3+v4)/2; miny=MAX(miny,-j);}
        if(v1 >= 0){ j=(v1+v4)/2; maxx=MIN(maxx,j);} }
  } else v4=-2;

// test for small enough area to fill

  if((maxx-minx) < CUTSV) {
    if(v2==1 && v4==0  ) { //A_b is diagonal
        FILLv2v4;
        continue;
    }
    if(v1==1 && v3==0 ) { //A_b is diagonal
        FILLv1v3;
        continue;
    }

    if(AREAbb < CUTS ) { // bb is small enough
        if(minx==0&&maxx==0){
            for(ly=miny; ly<=maxy; ly++) {edy=ly*ly; if( edy < dt[ly+y][x] ) dt[ly+y][x] = edy;}
        } else  {
            edy = minx*minx;
            for(ly=miny; ly<= maxy; ly++) {
                ed= edy+ly*ly;
                int *p = &dt[ly+y][minx+x];
                for(lx=minx; lx <=maxx; lx++,p++) { if(ed <*p) *p = ed; ed += 2*lx+1; }
            }
        }
        continue;
    }

  } // end of test on CUTSV

// search k{2} with k=1,2
  UNEAR12LV1; UNEAR12LV2;
//search k{1} with k=5 step D11
  DISTV12;

// test for small enough area to fill
  if(AREAbb < CUT12 ) {
    if(v2==1 && v4==0 ) {//A_b is diagonal
        FILLv2v4;
    } else if(v1==1 && v3==0  ) { //A_b is diagonal
        FILLv1v3;
    } else  {
        for(ly=miny; ly<= maxy; ly++) {
            VCUTF;
            OCUT(1);OCUT(2);
            int *p = &dt[ly+y][low+x]; ed = ly*ly + (low-1)*(low-1);
            for(lx=low; lx <=high; lx++,p++) { ed += 2*lx-1; if(ed <*p) *p = ed;}
        }
    }// end of filling
    continue;
  } // end of test on AREAbb < CUT12

// search k{4} with k=1
  UNEARV(3); UNEARV(4);
  UNEARV(5); UNEARV(6);
// search k{2} for k>2
  DISTOV(1); DISTOV(2);

// test for small enough area to fill
  if(AREAbb < CUT12L ) {
      if(v2==1 && v4==0 ) {
          FILLv2v4;
      } else if(v1==1 && v3==0  ) {
          FILLv1v3;
      } else {
          for(ly=miny; ly<= maxy; ly++) {
              VCUTF;
              OCUT(1);OCUT(2);
              int *p = &dt[ly+y][low+x]; ed = ly*ly + (low-1)*(low-1);
              for(lx=low; lx <=high; lx++,p++) { ed += 2*lx-1; if(ed <*p) *p = ed;}
          }
      }// end of filling
      continue;
  }

// search k{4} with k>1
  DISTOV(3); DISTOV(4); DISTOV(5); DISTOV(6);

// final filling
  if(v2==1 && v4==0) {
    FILLv2v4;
  } else if(v1==1 && v3==0) {
    FILLv1v3;
  } else {
    VCHECK;
    OCHECK(1);OCHECK(2);
    OCHECK(3);OCHECK(4); OCHECK(5);OCHECK(6);
    for(ly=miny; ly<= maxy; ly++) {//if(ly==0)continue; //leaving out is faster
        VCUTF;
        OCUT(1); OCUT(2);
        OCUT(3); OCUT(4); OCUT(5); OCUT(6);
        int *p = &dt[ly+y][low+x]; ed = ly*ly + (low-1)*(low-1); // incremental edy is slower
        for(lx=low; lx <=high; lx++,p++) { ed += 2*lx-1; if(ed <*p) *p = ed;}
    }
  }// end of filling


 } // end loop over bi
} // end of loop over y

std::cout<<"REUSSI"<<std::endl;
return 0;
}
