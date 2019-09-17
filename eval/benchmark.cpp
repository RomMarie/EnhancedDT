//****************
//* Romain MARIE *
//***************/

#include <ros/package.h>
#include <chrono>

// Différents critères
#include <distance_transform/static/eval/interface_eval_static_dt.h>
#include <distance_transform/static/eval/circles.h>
#include <distance_transform/static/eval/corners.h>
#include <distance_transform/static/eval/randomPixels.h>
#include <distance_transform/static/eval/randomSquares.h>
#include <distance_transform/static/eval/turningLine.h>
#include <distance_transform/static/eval/lenna.h>
#include <distance_transform/static/eval/kimia.h>

// Algos existants
#include <distance_transform/distance_transform.h>
int main(int argc, char ** argv){

    std::vector<dt::Algo*> algos; // Algos considérés
    std::vector<dt::eval::Critere*> criteres; // Critères considérés

    // On ajoute les différents algos
    algos.push_back(new dt::Felzenszwalb12());
    algos.push_back(new dt::Hirata96());
    algos.push_back(new dt::Lucet05());
    algos.push_back(new dt::Marie19());
    algos.push_back(new dt::Maurer03());
    algos.push_back(new dt::Schouten14());

    // On prépare le répertoire où sont stockés les résultats
    std::string path = ros::package::getPath("distance_transform");
    path += "/src/static/eval/";

    // On ajoute les différents critères
    criteres.push_back(new dt::eval::Corners(path,10,false));
    criteres.push_back(new dt::eval::Corners(path,10,true));
    criteres.push_back(new dt::eval::Circles(path,10,false));
    criteres.push_back(new dt::eval::Circles(path,10,true));
    criteres.push_back(new dt::eval::TurningLine(path,10,false));
    criteres.push_back(new dt::eval::TurningLine(path,10,true));
    criteres.push_back(new dt::eval::RandomPixels(path,10));
    criteres.push_back(new dt::eval::RandomSquares(path,10));
    criteres.push_back(new dt::eval::Lenna(path,10,false));
    criteres.push_back(new dt::eval::Lenna(path,10,true));
    criteres.push_back(new dt::eval::Kimia(path,10,false));
    criteres.push_back(new dt::eval::Kimia(path,10,true));
    criteres.push_back(new dt::eval::KimiaMean(path,10,false));
    criteres.push_back(new dt::eval::KimiaMean(path,10,true));

    // Pour chaque critère
    for(dt::eval::Critere* c:criteres){
        // Petit message pour indiquer l'état d'avancement du benchmark
        std::cout<<"\033[34mCritère : "<<c->getLabel()<<"\033[0m"<<std::endl;
        // Pour chaque algo
        for(dt::Algo* a:algos){
            // Petit message pour indiquer l'état d'avancement du benchmark
            std::cout<<"\033[33m"<<a->getLabel()<<"\033[0m"<<std::endl;
            // Exécution du test
            c->evaluer(a);
        }
    }

}
