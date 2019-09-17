//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/algo/maurer03.h>

dt::Maurer03::Maurer03(const cv::Mat &in):Algo("Maurer 2003")
{
    if(!in.empty()){
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);
        p_1d=cv::Mat(in.size(),cv::DataType<int>::type);
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
        g= std::vector<int>(static_cast<unsigned long>(in.cols));
        s= std::vector<int>(static_cast<unsigned long>(in.cols));
    }
}

void dt::Maurer03::prepare(const cv::Mat &in)
{
    // On prépare si besoin
    if(h.size()!=static_cast<unsigned long>(in.cols))
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(g.size()!=static_cast<unsigned long>(in.cols))
        g= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(s.size()!=static_cast<unsigned long>(in.cols))
        s= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(g_1d.rows!=in.rows||g_1d.cols!=in.cols)
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);
}

void dt::Maurer03::compute(const cv::Mat& in,cv::Mat& D,bool squared)
{
    // On prépare si besoin les matrices de sortie
    if(D.rows!=in.rows||D.cols!=in.cols)
        D=cv::Mat(in.size(),cv::DataType<double>::type);
    if(h.size()!=static_cast<unsigned long>(in.cols))
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(g.size()!=static_cast<unsigned long>(in.cols))
        g= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(s.size()!=static_cast<unsigned long>(in.cols))
        s= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(g_1d.rows!=in.rows||g_1d.cols!=in.cols)
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);

    // Puis on traite chaque dimension de manière séquentielle
    process1d(in);
    process2d(squared,D);
}

void dt::Maurer03::compute(const cv::Mat& in,cv::Mat& D,cv::Mat& P,bool squared)
{
    // On prépare si besoin les matrices de sortie
    if(D.rows!=in.rows||D.cols!=in.cols){
        D=cv::Mat(in.size(),cv::DataType<double>::type);
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
        g= std::vector<int>(static_cast<unsigned long>(in.cols));
        s= std::vector<int>(static_cast<unsigned long>(in.cols));
    }
    if(P.rows!=in.rows||P.cols!=in.cols)
        P=cv::Mat_<cv::Point>(in.size());
    if(g_1d.rows!=in.rows||g_1d.cols!=in.cols)
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);
    if(p_1d.rows!=in.rows||p_1d.cols!=in.cols)
        p_1d=cv::Mat(in.size(),cv::DataType<int>::type);

    // Puis on traite chaque dimension de manière séquentielle
    process1d_withP(in);
    process2d_withP(squared,D,P);
}

void dt::Maurer03::process1d(const cv::Mat &in)
{
    int rows=in.rows;
    int cols=in.cols;
    int maxD=rows+cols;

    for( int j=0;j<cols;j++){
        auto m = in.ptr<uchar>(0);
        auto g = g_1d.ptr<int>(0);
        if(m[j]==0)
            g[j]=0;
        else
            g[j]=maxD;
    }

    for(int i=1;i<rows;i++){
        auto m = in.ptr<uchar>(i);
        auto g = g_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(m[j]==0)
                g[j]=0;
            else
                g[j]=g_1d.at<int>(i-1,j)+1;
        }
    }

    for(int i=rows-2;i>=0;i--){
        auto g = g_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(g_1d.at<int>(i+1,j)<g[j]){
                g[j]=g_1d.at<int>(i+1,j)+1;
            }
        }
    }
}

void dt::Maurer03::process1d_withP(const cv::Mat &in)
{
    int rows=in.rows;
    int cols=in.cols;
    int maxD=rows+cols;

    for(int j=0;j<cols;j++){
        auto m = in.ptr<uchar>(0);
        auto g = g_1d.ptr<int>(0);
        if(m[j]==0)
            g[j]=0;
        else
            g[j]=maxD;
    }
    for(int i=1;i<rows;i++){
        auto m = in.ptr<uchar>(i);
        auto g = g_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(m[j]==0)
                g[j]=0;
            else
                g[j]=g_1d.at<int>(i-1,j)+1;
        }
    }

    for(int j=0;j<cols;j++){
        auto g = g_1d.ptr<int>(rows-1);
        auto p = p_1d.ptr<int>(rows-1);
        if(g[j]>=maxD)
            p[j]=-1;
        else
            p[j]=rows-1-g[j];
    }

    for(int i=rows-2;i>=0;i--){
        auto g = g_1d.ptr<int>(i);
        auto p = p_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(g_1d.at<int>(i+1,j)<g[j]){
                g[j]=g_1d.at<int>(i+1,j)+1;
                p[j]=g[j]+i;
            }
            else if(g[j]>=maxD)
                p[j]=-1;
            else
                p[j]=i-g[j];
        }
    }
}

void dt::Maurer03::process2d(bool squared, cv::Mat &D)
{
    int rows=D.rows;
    int cols=D.cols;
    int maxD=rows+cols;
    int i,j;
    int l, f1, ns, t0;

#define REM(glm1,gl,f,hlm1,hl,x) ( (x-hlm1)*(gl - (x-hl)*(hl-hlm1) ) -(x-hl)*glm1-(hl-hlm1)*f  >0)

    for(i=0; i < rows; i++) {
        auto pG=g_1d.ptr<int>(i);
        l=-1;

        // Optimisation 1
        for(j=0;j<cols;j++) {
            s[static_cast<unsigned long>(j)]=pG[j]*pG[j];
        }
        for(j=0;j<cols;j++) {
            f1=s[static_cast<unsigned long>(j)];
            // Optimisation 2
            if(pG[j]>=maxD)
                continue;


            if(l<1){
                l++; g[static_cast<unsigned long>(l)]=f1; h[static_cast<unsigned long>(l)]= j;
            }
            else {
                while (l >=1 && REM(g[static_cast<unsigned long>(l)-1],g[static_cast<unsigned long>(l)],f1,h[static_cast<unsigned long>(l)-1],h[static_cast<unsigned long>(l)],j) )
                    l--;
                l++; g[static_cast<unsigned long>(l)]=f1; h[static_cast<unsigned long>(l)]=j;
            }

        }

        auto pD = D.ptr<double>(i);
        if(l==0) {
            for(j=0;j<cols;j++)
                if(squared)
                    pD[j]=g[static_cast<unsigned long>(l)]+(h[static_cast<unsigned long>(l)]-j)*(h[static_cast<unsigned long>(l)]-j);
                else
                    pD[j]=sqrt(g[static_cast<unsigned long>(l)]+(h[static_cast<unsigned long>(l)]-j)*(h[static_cast<unsigned long>(l)]-j));
        }
        else if(l>-1) {
            ns=l; l=0;
            for(j=0;j<cols;j++) {
                t0= g[static_cast<unsigned long>(l)]+(h[static_cast<unsigned long>(l)]-j)*(h[static_cast<unsigned long>(l)]-j);
                while (l < ns && ( t0 > g[static_cast<unsigned long>(l)+1]+(h[static_cast<unsigned long>(l)+1]-j)*(h[static_cast<unsigned long>(l)+1]-j) ) ) {
                    t0=g[static_cast<unsigned long>(l)+1]+(h[static_cast<unsigned long>(l)+1]-j)*(h[static_cast<unsigned long>(l)+1]-j); l++;
                }
                if(squared)
                    pD[j]= t0;
                else
                    pD[j]=sqrt(t0);
            }
        }
        else{
            for(j=0;j<cols;j++) {
                if(squared)
                    pD[j]= maxD*maxD;
                else
                    pD[j]=maxD;
            }
        }
    }
}

void dt::Maurer03::process2d_withP(bool squared, cv::Mat &D, cv::Mat &P)
{
    int rows=D.rows;
    int cols=D.cols;
    int maxD=rows+cols;
    int i,j;
    int l, f1, ns, t0;

#define REM(glm1,gl,f,hlm1,hl,x) ( (x-hlm1)*(gl - (x-hl)*(hl-hlm1) ) -(x-hl)*glm1-(hl-hlm1)*f  >0)

    for(i=0; i < rows; i++) {
        auto pG=g_1d.ptr<int>(i);
        l=-1;

        // Optimisation 1
        for(j=0;j<cols;j++) {
            s[static_cast<unsigned long>(j)]=pG[j]*pG[j];
        }
        for(j=0;j<cols;j++) {
            f1=s[static_cast<unsigned long>(j)];
            // Optimisation 2
            if(pG[j]>=maxD)
                continue;

            if(l<1){
                l++; g[static_cast<unsigned long>(l)]=f1; h[static_cast<unsigned long>(l)]= j;
            }
            else {
                while (l >=1 && REM(g[static_cast<unsigned long>(l)-1],g[static_cast<unsigned long>(l)],f1,h[static_cast<unsigned long>(l)-1],h[static_cast<unsigned long>(l)],j) )
                    l--;
                l++; g[static_cast<unsigned long>(l)]=f1;h[static_cast<unsigned long>(l)]=j;
            }

        }
        auto pD = D.ptr<double>(i);
        auto pP = P.ptr<cv::Point>(i);
        auto pP1 = p_1d.ptr<int>(i);
        if(l==0) {
            for(j=0;j<cols;j++){
                if(squared)
                    pD[j]=g[static_cast<unsigned long>(l)]+(h[static_cast<unsigned long>(l)]-j)*(h[static_cast<unsigned long>(l)]-j);
                else
                    pD[j]=sqrt(g[static_cast<unsigned long>(l)]+(h[static_cast<unsigned long>(l)]-j)*(h[static_cast<unsigned long>(l)]-j));
                pP[j]=cv::Point(pP1[h[static_cast<unsigned long>(l)]],h[static_cast<unsigned long>(l)]);
            }
        }
        else if(l>-1) {
            ns=l; l=0;
            for(j=0;j<cols;j++) {
                t0= g[static_cast<unsigned long>(l)]+(h[static_cast<unsigned long>(l)]-j)*(h[static_cast<unsigned long>(l)]-j);
                while (l < ns && ( t0 > g[static_cast<unsigned long>(l)+1]+(h[static_cast<unsigned long>(l)+1]-j)*(h[static_cast<unsigned long>(l)+1]-j) ) ) {
                    t0=g[static_cast<unsigned long>(l)+1]+(h[static_cast<unsigned long>(l)+1]-j)*(h[static_cast<unsigned long>(l)+1]-j); l++;
                }
                if(squared)
                    pD[j]= t0;
                else
                    pD[j]=sqrt(t0);
                pP[j]=cv::Point(pP1[h[static_cast<unsigned long>(l)]],h[static_cast<unsigned long>(l)]);
            }
        }
        else{
            for(j=0;j<cols;j++) {
                if(squared)
                    pD[j]= maxD*maxD;
                else
                    pD[j]=maxD;
                pP[j]=cv::Point(-1,-1);
            }
        }
    }
}

