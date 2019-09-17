//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/eval/turningLine.h>
#include <chrono>
dt::eval::TurningLine::TurningLine(std::string path, int iter, bool inverse):
    Critere("TurningLine"+std::string(inverse?"Inverse":""),path+"resultats/",iter),img(36)
{
    fichier<<"-- Turning Line --";
    for(unsigned long j=0;j<36;j++){
        std::ostringstream ss;
        ss<<path<<"images/turning_line/s1000/";
        ss<<"i"<<std::setfill('0')<<std::setw(3)<<j<<".png";
        img[j]=cv::imread(ss.str(),0);
        if(inverse)
            img[j]=cv::Scalar::all(255)-img[j];
    }

    for(int j=0;j<36;j++){
        fichier<<","<<j*5<<"°";
    }
    fichier<<std::endl;
}

void dt::eval::TurningLine::evaluer(dt::Algo *algo)
{
    double step = 100/36.;
    double pourcent = 0;
    fichier<<algo->getLabel()<<",";
    for(unsigned long j=0;j<36;j++){ // Pour chaque taille
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
