//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/eval/randomSquares.h>
#include <chrono>
dt::eval::RandomSquares::RandomSquares(std::string path, int iter, bool inverse):
    Critere("RandomSquares"+std::string(inverse?"Inverse":""),path+"resultats/",iter),img(7)
{
    fichier<<"-- Random Squares --";
    for(unsigned long i=0;i<7;i++){
        img[i]=std::vector<cv::Mat>(5);
        for(unsigned long j=0;j<5;j++){

            std::ostringstream ss;
            ss<<path<<"images/random_squares/3000/"<<i*15<<"/";
            ss<<"i"<<std::setfill('0')<<std::setw(3)<<j<<".png";
            img[i][j]=cv::imread(ss.str(),0);
        }
    }

    fichier<<",10%,25%,50%,75%,90%";
    fichier<<std::endl;
}

void dt::eval::RandomSquares::evaluer(dt::Algo *algo)
{
    double step = 100/35.;
    double pourcent = 0;
    fichier<<algo->getLabel();
    for(unsigned long i=0;i<7;i++){
        for(unsigned long j=0;j<5;j++){ // Pour chaque taille
            cv::Mat D=cv::Mat::zeros(img[i][j].size(),cv::DataType<double>::type);
            algo->prepare(img[i][j]);

            std::chrono::steady_clock::time_point _start(std::chrono::steady_clock::now());
            for(int k=0;k<iter;k++){ // Pour chaque iteration
                algo->compute(img[i][j],D,true);
            }
            std::chrono::steady_clock::time_point _end(std::chrono::steady_clock::now());
            pourcent+=step;
            std::cout<<'\r'<<static_cast<int>(pourcent)<<"/100"<<std::flush;
            fichier<<","<<std::chrono::duration_cast<std::chrono::duration<double>>(_end - _start).count()/(iter*D.rows*D.cols)*1000000000;
        }
        fichier<<std::endl;
    }

    std::cout<<'\r';
}
