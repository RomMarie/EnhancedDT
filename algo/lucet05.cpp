//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/algo/lucet05.h>

dt::Lucet05::Lucet05(const cv::Mat &in):Algo("Lucet 2005")
{
}



int lowhull(const int n, const int * x, const int* y, int *xh, int *yh){
    int v; // index of last point on the hull
    int i;
    int a,b; // for comparing an input point with the last two points on the hull, two parts to avoid division

    xh[0]=x[0]; yh[0] = y[0]; // copy first input point
    if(n==1) return 1;
    xh[1]=x[1]; yh[1] = y[1]; // copy second input point
    v=1;
    for(i=2; i < n; i++) { // loop over the rest of the input point
        a= (y[i]-yh[v])*(xh[v]-xh[v-1]); b= (yh[v]-yh[v-1]) * (x[i]-xh[v]);
        // optimization, a few percent faster:
        if( a > b) { // often occuring case
            xh[v+1]=x[i]; yh[v+1]=y[i]; v++; continue; // add point to hull
        }
        if(v>1) {  v--; a= (y[i]-yh[v])*(xh[v]-xh[v-1]); b= (yh[v]-yh[v-1]) * (x[i]-xh[v]);} // remove last point from hull
        else { xh[1]=x[i]; yh[1]=y[i]; continue;} // only 2 points on hull
        //end of optimization
        for(;v>1;) { // a while loop while(v>1 && a<=b) is slower
            if(a>b) {  break;} // add point to hull
            v--; a= (y[i]-yh[v])*(xh[v]-xh[v-1]); b= (yh[v]-yh[v-1]) * (x[i]-xh[v]); // remove last point from hull
        }
        if( a > b  ) { xh[v+1]=x[i]; yh[v+1]=y[i];v++; }
        else { xh[1]=x[i]; yh[1]=y[i]; } // only 2 points on hull

    }
    return v+1;
}

void conjpart(const int abl, const int * ab, const int* abval, int* result, const int length) {

    int i,j;

    // specific case if only 1 point on the hull
    if (abl==1){
        for(i=0; i<length;i++) result[i]=(i+1)*( 2*(ab[0]+1))-abval[0]; // factor 2 to compensate not doing the division by 2 for abval[]
        return ;
    } // otherwise we have at least 2 points on the convex hull, and so at least one slope

    int a,b;
    i=1; // first slope
    a = abval[1] - abval[0] ; b = ab[1]-ab[0] ; // first slope c_0=a/b
    for(j=1;j<=length;){
        if ((2*j)*b <=a) { result[j-1]=((2*j)*(ab[i-1]+1))-abval[i-1]; j++; } // often occuring case, note factors 2
        else {
            if(i== abl-1){for(;j<=length;j++) result[j-1]=((2*j)*(ab[i]+1))-abval[i];} // directly go to end of j
            else {a = abval[i+1] - abval[i] ; b = ab[i+1]-ab[i] ; i++;} // next slope
        }
    }

    return ;
}

void dt::Lucet05::compute(const cv::Mat& in,cv::Mat& D,bool squared)
{
    if(D.rows!=in.rows||D.cols!=in.cols){
        D=cv::Mat(in.size(),cv::DataType<double>::type);
    }

    int rows = in.rows,cols = in.cols, maxD=(rows*cols)*(rows*cols);

    int i, j, k;  //loop indices

    if(rows > cols)
        k=rows;
    else
        k=cols;
    int *x = new int [k];
    int *y = new int [k];
    int *xh = new int [k];
    int *yh = new int [k];
    int *result = new int [k]; // input and result vectors
    int *result2 = new int [k]; // input and result vectors

    // first dimension

    for(i=0;i<rows;i++){
        auto m = in.ptr<uchar>(i);
        auto d = D.ptr<double>(i);
        k=0; // number of object pixels in row
        for(j=0;j<cols;j++){
            if (m[j]==0){ // object pixel
                y[k]= (i+1)*(i+1)+(j+1)*(j+1); // do not use division by 2, compensate this in conjpart()
                x[k]=j; k++;
            }
        }
        if(k==0){
            for(j=0;j<cols;j++) {
                d[j] = maxD;
            }
        }
        else{
            //std::cout<<"ici"<<" ";
            conjpart(k,x, y,result2, cols);  // lowhull() is not needed in first dimension, output of conjpart() stored in output matrix

            for(j=0;j<cols;j++) {
                d[j] = result2[j];
            }
        }
    }

    // second dimension
    for(j=0;j<cols;j++){
        k=0;
        for(i=0;i<rows;i++){
            auto d = D.ptr<double>(i);
            if( d[j] < maxD) {
                y[k]= -static_cast<int>(d[j]); x[k]=i; k++;
            }
        }

        i= lowhull(k,x,y,xh,yh);
        conjpart(i,xh,yh, result, rows);
        for(i=0;i<rows;i++){
            auto d = D.ptr<double>(i);
            if(result[i]>(i+1)*(i+1)+ (j+1)*(j+1))// if NaN
                d[j]=maxD;
            else
                d[j] = (i+1)*(i+1)+ (j+1)*(j+1) - result[i] ;
        } // note  that factor 2 is not in here
    }

    delete []x; delete []y;
    delete []xh; delete []yh; delete []result; delete []result2;

    if(!squared){
        for(i=0;i<rows;i++){
            auto d = D.ptr<double>(i);
            for(j=0;j<cols;j++){
                d[j]=sqrt(d[j]);
            }
        }
    }
}


// Non implémenté !!!!!
void dt::Lucet05::compute(const cv::Mat& in,cv::Mat& D,cv::Mat& P,bool squared)
{

    if(D.rows!=in.rows||D.cols!=in.cols){
        D=cv::Mat(in.size(),cv::DataType<double>::type);
    }
    if(P.rows!=in.rows||P.cols!=in.cols){
        P=cv::Mat_<cv::Point>::zeros(in.size());
    }

    compute(in,D,squared);

}
