//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/eval/lenna.h>
#include <chrono>
dt::eval::Lenna::Lenna(std::__cxx11::string path, int iter, bool inverse):
    Critere("Lenna"+std::string(inverse?"Inverse":""),path+"resultats/",iter),img(4)
{
    fichier<<"-- Lenna --";
    for(unsigned long j=0;j<4;j++){
        std::ostringstream ss;
        ss<<path<<"images/Lenna_edges.png";
        img[j]=cv::imread(ss.str(),0);
        if(inverse)
            img[j]=cv::Scalar::all(255)-img[j];
        cv::resize(img[j],img[j],img[j].size()*static_cast<int>(j+1),0,0,cv::INTER_NEAREST);
    }

    for(unsigned long j=0;j<4;j++){
        fichier<<","<<img[j].rows*img[j].cols<<"%";
    }
    fichier<<std::endl;
}

void dt::eval::Lenna::evaluer(dt::Algo *algo)
{
    double step = 25;
    double pourcent = 0;
    fichier<<algo->getLabel()<<",";
    for(unsigned long j=0;j<4;j++){ // Pour chaque taille
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
