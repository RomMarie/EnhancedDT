//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/eval/randomPixels.h>
#include <chrono>
dt::eval::RandomPixels::RandomPixels(std::string path, int iter):
    Critere("RandomPixels",path+"resultats/",iter),img(11)
{
    fichier<<"-- Random Pixels --";
    for(unsigned long j=0;j<11;j++){
        std::ostringstream ss;
        ss<<path<<"images/random_points/";
        ss<<"i"<<std::setfill('0')<<std::setw(3)<<j<<".png";
        img[j]=cv::imread(ss.str(),0);
    }

    for(int j=0;j<11;j++){
        fichier<<","<<j*9<<"%";
    }
    fichier<<std::endl;
}

void dt::eval::RandomPixels::evaluer(dt::Algo *algo)
{
    double step = 100/11.;
    double pourcent = 0;
    fichier<<algo->getLabel()<<",";
    for(unsigned long j=0;j<11;j++){ // Pour chaque taille
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
