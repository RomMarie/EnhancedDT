//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/algo/hirata96.h>


#define F(y,j) ( h[static_cast<unsigned long>(y)] + ((j)-(y)) * ((j)-(y)) )
#define SEP(u,v) ( (v)*(v) -(u)*(u) + h[static_cast<unsigned long>(v)] - h[static_cast<unsigned long>(u)]) / ( 2.* ( (v)-(u) ) )

dt::Hirata96::Hirata96(const cv::Mat &in):Algo("Hirata 1996")
{
    if(!in.empty()){
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);
        p_1d=cv::Mat(in.size(),cv::DataType<int>::type);
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
        obs=std::vector<obst>(static_cast<unsigned long>(in.cols)+1);
    }
}

void dt::Hirata96::prepare(const cv::Mat &in)
{
    // On prépare si besoin
    if(h.size()!=static_cast<unsigned long>(in.cols))
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(obs.size()!=static_cast<unsigned long>(in.cols)+1)
        obs= std::vector<obst>(static_cast<unsigned long>(in.cols)+1);
    if(g_1d.rows!=in.rows||g_1d.cols!=in.cols)
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);
}
void dt::Hirata96::compute(const cv::Mat& in,cv::Mat& D,bool squared)
{
    // On prépare si besoin les matrices de sortie
    if(D.rows!=in.rows||D.cols!=in.cols)
        D=cv::Mat(in.size(),cv::DataType<double>::type);
    if(h.size()!=static_cast<unsigned long>(in.cols))
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(obs.size()!=static_cast<unsigned long>(in.cols)+1)
        obs=std::vector<obst>(static_cast<unsigned long>(in.cols)+1);
    if(g_1d.rows!=in.rows||g_1d.cols!=in.cols)
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);

    // Puis on traite chaque dimension de manière séquentielle
    process1d(in);
    process2d(squared,D);
}

void dt::Hirata96::compute(const cv::Mat& in,cv::Mat& D,cv::Mat& P,bool squared)
{
    // On prépare si besoin les matrices de sortie
    if(D.rows!=in.rows||D.cols!=in.cols)
        D=cv::Mat(in.size(),cv::DataType<double>::type);
    if(h.size()!=static_cast<unsigned long>(in.cols))
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(obs.size()!=static_cast<unsigned long>(in.cols)+1)
        obs=std::vector<obst>(static_cast<unsigned long>(in.cols)+1);
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

void dt::Hirata96::process1d(const cv::Mat &in)
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

    for(int i=rows-2;i>=0;i--){
        auto g = g_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(g_1d.at<int>(i+1,j)<g[j]){
                g[j]=g_1d.at<int>(i+1,j)+1;
            }
        }
    }
}

void dt::Hirata96::process1d_withP(const cv::Mat &in)
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

void dt::Hirata96::process2d(bool squared, cv::Mat &D)
{

    int rows=D.rows, cols=D.cols, maxD=rows+cols;
    int w;
    for(int i=0; i < rows; i++) {
        auto q=obs.begin()+1;
        q->s=0; q->t=0;
        auto pG=g_1d.ptr<int>(i);

        // Optimisation 1
        for(int j=0; j<cols;j++)
            h[static_cast<unsigned long>(j)]= pG[j]*pG[j];

        for(int j=1;j<cols;j++) {
            // Optimisation 2
            if(pG[j]>=maxD)
                continue;

            while ( q !=obs.begin() and F(q->s,q->t) > F(j,q->t) )
                q--;
            if( q == obs.begin()){
                q++; q->s=j;
            }
            else {
                w = static_cast<int>(1+ SEP(q->s,j));
                if( w < cols ){
                    q++; q->s=j; q->t= w;
                }
            }
        }

        auto pD = D.ptr<double>(i);
        for(int j=cols-1; j>=0; j--) {
            if(squared)
                pD[j] = F(q->s,j);
            else
                pD[j] = sqrt(F(q->s,j));

            if( j == q->t )
                q--;
        }
    }
}

void dt::Hirata96::process2d_withP(bool squared, cv::Mat &D, cv::Mat &P)
{

    int rows=D.rows, cols=D.cols, maxD=rows+cols;
    int w;
    for(int i=0; i < rows; i++) {
        auto q=obs.begin()+1;
        q->s=0; q->t=0;
        auto pG=g_1d.ptr<int>(i);

        // Optimisation 1
        for(int j=0; j<cols;j++)
            h[static_cast<unsigned long>(j)]= pG[j]*pG[j];

        for(int j=1;j<cols;j++) {
            // Optimisation 2
            if(pG[j]>=maxD)
                continue;

            while ( q != obs.begin() and F(q->s,q->t) > F(j,q->t) )
                q--;
            if( q == obs.begin()){
                q++; q->s=j;
            }
            else {
                w = static_cast<int>(1+ SEP(q->s,j));
                if( w < cols ){
                    q++; q->s=j; q->t= w;
                }
            }
        }

        auto pD = D.ptr<double>(i);
        auto pP = P.ptr<cv::Point>(i);
        auto pP1 = p_1d.ptr<int>(i);
        for(int j=cols-1; j>=0; j--) {
            if(squared)
                pD[j] = F(q->s,j);
            else
                pD[j] = sqrt(F(q->s,j));
            pP[j]=cv::Point(pP1[q->s],q->s);
            if( j == q->t )
                q--;
        }
    }
}

