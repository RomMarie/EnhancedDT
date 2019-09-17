//****************
//* Romain MARIE *
//***************/

// Juste pour comparaison dans article 2019 ... L'algo semble ne pas marcher dans
// tous les cas de figure, ne colle pas avec l'interface, et est TRES dur à configurer
// Ne pas utiliser dans application concrète !

#include <distance_transform/static/algo/schouten14.h>

dt::Schouten14::Schouten14(const cv::Mat &in):Algo("Schouten 2014"),
    feed(in.rows,in.cols),
    inCImg(static_cast<unsigned int>(in.cols),static_cast<unsigned int>(in.rows)),
    outCImg(static_cast<unsigned int>(in.cols),static_cast<unsigned int>(in.rows))
{

}

void dt::Schouten14::prepare(const cv::Mat &in)
{
    if(inCImg.height!=static_cast<unsigned int>(in.rows)||inCImg.width!=static_cast<unsigned int>(in.cols)){
        inCImg = inCImg.assign(static_cast<unsigned int>(in.cols),static_cast<unsigned int>(in.rows));
        outCImg = outCImg.assign(static_cast<unsigned int>(in.cols),static_cast<unsigned int>(in.rows));
    }
    uchar* d = inCImg.data;
    for(unsigned int i=0;i<inCImg.height;i++){
        for(unsigned int j=0;j<inCImg.width;j++){
            if(i<4||j<4||i>=inCImg.height-4||j>=inCImg.width-4)
                *(d++)=255;
            else
                *(d++)=in.at<uchar>(static_cast<int>(i),static_cast<int>(j));

        }
    }

}

void dt::Schouten14::finalize(const cv::Mat& in,cv::Mat &D)
{
    D = cv::Mat::zeros(in.size(),cv::DataType<double>::type);
    for(int i=0;i<D.rows;i++){
        for(int j=0;j<D.cols;j++){
                D.at<double>(i,j)=outCImg.at(static_cast<unsigned int>(j),static_cast<unsigned int>(i));
        }
    }
}

void dt::Schouten14::compute(const cv::Mat& in,cv::Mat& D,bool squared)
{
    feed.genfeed(inCImg,outCImg);
    if(!squared)
        outCImg=outCImg.sqrt();

}

void dt::Schouten14::compute(const cv::Mat& in,cv::Mat& D,cv::Mat& P,bool squared)
{
}
