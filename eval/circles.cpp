//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/eval/circles.h>
#include <chrono>
dt::eval::Circles::Circles(std::string path, int iter, bool inverse):
    Critere("Cercles"+std::string(inverse?"Inverse":""),path+"resultats/",iter),img(8)
{

    fichier<<"-- Circles --";
    // 8 tailles de disques sont considérées
    for(unsigned long j=0;j<8;j++){

        // On récupère et on stocke l'image correspondant à la taille courante
        std::stringstream ss;
        ss<<path<<"images/inscribed_circle/";
        ss<<"i00"<<j<<".png";
        img[j]=cv::imread(ss.str(),0);

        // Et on l'inverse si demandée
        if(inverse)
            img[j]=cv::Scalar::all(255)-img[j];
    }

    // Pour chaque taille, on écrit l'intitulé de la colonne résultat correspondant
    for(unsigned long j=0;j<8;j++){
        fichier<<","<<img[j].rows*img[j].rows;
    }
    fichier<<std::endl;
}

void dt::eval::Circles::evaluer(dt::Algo *algo)
{
    double step = 100/8.;
    double pourcent = 0; // Variable pour afficher l'état courant d'exécution du critère

    // La ligne courante du fichier résultat correspond à un algo
    fichier<<algo->getLabel()<<",";
    for(unsigned long j=0;j<8;j++){ // Pour chaque taille
        cv::Mat D=cv::Mat::zeros(img[j].size(),cv::DataType<double>::type);
        algo->prepare(img[j]); // Certains algos en ont besoin

        // On démarre le chrono
        std::chrono::steady_clock::time_point _start(std::chrono::steady_clock::now());
        for(int k=0;k<iter;k++){ // Pour chaque iteration
            // Ici, seuls les calculs sont mesurés, pas les initialisations mémoire
            algo->compute(img[j],D,true);
        }
        // On arrête le chrono après l'ensemble des itérations
        std::chrono::steady_clock::time_point _end(std::chrono::steady_clock::now());
        pourcent+=step;
        std::cout<<'\r'<<static_cast<int>(pourcent)<<"/100"<<std::flush;

        // Le temps retenu est la moyenne par pixel (exprimé en ns)
        fichier<<std::chrono::duration_cast<std::chrono::duration<double>>(_end - _start).count()/(iter*D.rows*D.cols)*1000000000<<",";
    }
    // On passe à la ligne dans l'attente du prochain algo
    fichier<<std::endl;

    std::cout<<'\r';
}
