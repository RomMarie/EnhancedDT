//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/eval/kimia.h>
#include <chrono>
dt::eval::Kimia::Kimia(std::string path, int iter, bool inverse):
    Critere("Kimia"+std::string(inverse?"Inverse":""),path+"resultats/",iter),img(216)
{
    fichier<<"-- Kimia --";
    for(unsigned long j=1;j<217;j++){
        std::ostringstream ss;
        ss<<path<<"images/kimia/";
        ss<<"img ("<<j<<").pgm";
        img[j-1]=cv::imread(ss.str(),0);
        //cv::resize(img[j-1],img[j-1],cv::Size(img[j-1].cols*4,img[j-1].rows*4),0,0,cv::INTER_NEAREST);

        if(!inverse)
            img[j-1]=cv::Scalar::all(255)-img[j-1];
    }

    for(int j=0;j<216;j++){
        fichier<<","<<j+1;
    }
    fichier<<std::endl;
}

void dt::eval::Kimia::evaluer(dt::Algo *algo)
{
    double step = 100/216.;
    double pourcent = 0;
    fichier<<algo->getLabel()<<",";
    for(unsigned long j=0;j<216;j++){ // Pour chaque taille
        cv::Mat D=cv::Mat::zeros(img[j].size(),cv::DataType<double>::type);
        algo->prepare(img[j]);
        std::chrono::steady_clock::time_point _start(std::chrono::steady_clock::now());
        for(int k=0;k<iter;k++){ // Pour chaque iteration
            algo->compute(img[j],D,true);
        }
        std::chrono::steady_clock::time_point _end(std::chrono::steady_clock::now());
        pourcent+=step;
        std::cout<<'\r'<<static_cast<int>(pourcent)<<"/100"<<std::flush;
        fichier<<std::chrono::duration_cast<std::chrono::duration<double>>(_end - _start).count()/(iter*D.rows*D.cols)*1000000000<<",";
    }
    fichier<<std::endl;

    std::cout<<'\r';
}

dt::eval::KimiaMean::KimiaMean(std::string path, int iter, bool inverse):
    Critere("KimiaMean"+std::string(inverse?"Inverse":""),path+"resultats/",iter),img(216)
{
    fichier<<"-- Kimia Mean--";
    for(unsigned long j=1;j<217;j++){
        std::ostringstream ss;
        ss<<path<<"images/kimia/";
        ss<<"img ("<<j<<").pgm";
        img[j-1]=cv::imread(ss.str(),0);
        //cv::resize(img[j-1],img[j-1],cv::Size(img[j-1].cols*4,img[j-1].rows*4),0,0,cv::INTER_NEAREST);

        if(!inverse)
            img[j-1]=cv::Scalar::all(255)-img[j-1];
    }

    for(int j=0;j<10;j++){
        fichier<<","<<j+1;
    }
    fichier<<std::endl;
}

void dt::eval::KimiaMean::evaluer(dt::Algo *algo)
{
    double step = 10.;
    double pourcent = 0;
    fichier<<algo->getLabel()<<",";

    for(int t = 1;t<=10;t++){
        double temps = 0;
        for(unsigned long j=0;j<216;j++){ // Pour chaque taille
            cv::Mat im;
            cv::resize(img[j],im,cv::Size(img[j].cols*t,img[j].rows*t),0,0,cv::INTER_NEAREST);

            cv::Mat D=cv::Mat::zeros(im.size(),cv::DataType<double>::type);
            algo->prepare(im);
            std::chrono::steady_clock::time_point _start(std::chrono::steady_clock::now());
            for(int k=0;k<iter;k++){ // Pour chaque iteration
                algo->compute(im,D,true);
            }
            std::chrono::steady_clock::time_point _end(std::chrono::steady_clock::now());
            temps+=std::chrono::duration_cast<std::chrono::duration<double>>(_end - _start).count()/(iter*D.rows*D.cols)*1000000000;
        }
        pourcent+=step;
        std::cout<<'\r'<<static_cast<int>(pourcent)<<"/100"<<std::flush;
        fichier<<temps/216.<<",";
    }
    fichier<<std::endl;
    std::cout<<'\r'<<"";
}
