//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/eval/corners.h>
#include <chrono>
dt::eval::Corners::Corners(std::__cxx11::string path, int iter, bool inverse):
    Critere("Corners"+std::string(inverse?"Inverse":""),path+"resultats/",iter),img(4)
{
    fichier<<"-- Corners --";
    // 4 coins
    for(unsigned long i=0;i<4;i++){
        img[i]=std::vector<cv::Mat>(10);
        // 10 tailles différentes
        for(unsigned long j=0;j<10;j++){

            // On construit les images
            if(inverse){ // pixel blanc sur fond noir
                img[i][j]=cv::Mat::zeros(static_cast<int>(j+1)*200,static_cast<int>(j+1)*200,cv::DataType<uchar>::type);

                // 4,4 car Lucet05 et Schouten14 en ont besoin
                img[i][j].at<uchar>(4,4)=255;
                // Rectangle contour car Lucet05 et Schouten14 en ont besoin
                cv::rectangle(img[i][j],cv::Point(0,0),cv::Point(static_cast<int>(j+1)*200-1,static_cast<int>(j+1)*200-1),255,4);
            }
            else{ // pixel noir sur fond blanc
                img[i][j]=cv::Mat::ones(static_cast<int>(j+1)*200,static_cast<int>(j+1)*200,cv::DataType<uchar>::type)*255;
                // 4,4 car Lucet05 et Schouten14 en ont besoin
                img[i][j].at<uchar>(4,4)=0;
                // Rectangle contour car Lucet05 et Schouten14 en ont besoin
                cv::rectangle(img[i][j],cv::Point(0,0),cv::Point(static_cast<int>(j+1)*200-1,static_cast<int>(j+1)*200-1),0,4);
            }
        }
    }

    // Pour chaque taille
    for(unsigned long j=0;j<10;j++){
        // On donne l'intitulé de la colonne résultat correspondant
        fichier<<","<<img[0][j].rows*img[0][j].rows;
    }
    fichier<<std::endl;
}

void dt::eval::Corners::evaluer(dt::Algo *algo)
{
    // On anote la ligne dans le fichier résultat
    fichier<<algo->getLabel()<<std::endl;

    double step = 2.5;
    double pourcent = 0;
    for(unsigned long i=0;i<4;i++){ // Pour chaque coin
        fichier<<algo->getLabel()<<": Coin"<<i+1<<",";
        for(unsigned long j=0;j<10;j++){ // Pour chaque taille
            cv::Mat D=cv::Mat::zeros(img[i][j].size(),cv::DataType<double>::type);
            algo->prepare(img[i][j]);

            std::chrono::steady_clock::time_point _start(std::chrono::steady_clock::now());
            for(int k=0;k<iter;k++){ // Pour chaque iteration
                algo->compute(img[i][j],D,true);
            }
            std::chrono::steady_clock::time_point _end(std::chrono::steady_clock::now());

            pourcent+=step;
            std::cout<<'\r'<<static_cast<int>(pourcent)<<"/100"<<std::flush;
            fichier<<std::chrono::duration_cast<std::chrono::duration<double>>(_end - _start).count()/(iter*D.rows*D.cols)*1000000000<<",";
        }
        fichier<<std::endl;
    }
    std::cout<<'\r';
}
