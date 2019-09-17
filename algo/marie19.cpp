//****************
//* Romain MARIE *
//***************/

#include <distance_transform/static/algo/marie19.h>
#include <ros/ros.h>
#include <thread>

// Fonctions de l'algo. Préprocesseur pour efficacité !
#define F(y,j) ( h[static_cast<unsigned long>(y)] + ((j)-(y)) * ((j)-(y)) )
#define SEP(u,v) ( (v)*(v) -(u)*(u) + h[static_cast<unsigned long>(v)] - h[static_cast<unsigned long>(u)]) / ( 2.* ( (v)-(u) ) )
#define SEP2(u,v) (v-sqrt(u*u-h[static_cast<unsigned long>(v)]))
#define SEP3(u,v) (v-sqrt(u-h[static_cast<unsigned long>(v)])) // Comme SEP2, mais u est passé au carré


dt::Marie19::Marie19(const cv::Mat &in):Algo("Marie 2019"){
    if(!in.empty()){ // Si l'image témoin est fournie, on alloue la mémoire
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);
        p_1d=cv::Mat(in.size(),cv::DataType<int>::type);

        h= std::vector<int>(static_cast<unsigned long>(in.cols));
        obs=std::vector<obst>(static_cast<unsigned long>(in.cols)+1);
    }
}

void dt::Marie19::prepare(const cv::Mat &in)
{
    // On prépare si besoin
    if(h.size()!=static_cast<unsigned long>(in.cols))
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(obs.size()!=static_cast<unsigned long>(in.cols)+1)
        obs= std::vector<obst>(static_cast<unsigned long>(in.cols)+1);
    if(g_1d.rows!=in.rows||g_1d.cols!=in.cols)
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);
    if(p_1d.rows!=in.rows||p_1d.cols!=in.cols)
        p_1d=cv::Mat(in.size(),cv::DataType<int>::type);
}

void dt::Marie19::compute(const cv::Mat &in, cv::Mat &D, bool squared)
{
    // On prépare si besoin les matrices de sortie
    if(D.rows!=in.rows||D.cols!=in.cols)
        D=cv::Mat(in.size(),cv::DataType<double>::type);
    if(h.size()!=static_cast<unsigned long>(in.cols))
        h= std::vector<int>(static_cast<unsigned long>(in.cols));
    if(obs.size()!=static_cast<unsigned long>(in.cols)+1)
        obs= std::vector<obst>(static_cast<unsigned long>(in.cols)+1);
    if(g_1d.rows!=in.rows||g_1d.cols!=in.cols)
        g_1d=cv::Mat(in.size(),cv::DataType<int>::type);

    // Puis on traite chaque dimension de manière séquentielle
    process1d(in);
    process2d(squared,D);
}

void dt::Marie19::compute(const cv::Mat &in, cv::Mat &D, cv::Mat &P, bool squared)
{
    // On prépare si besoin les matrices de sortie
    if(D.rows!=in.rows||D.cols!=in.cols){
        D=cv::Mat(in.size(),cv::DataType<double>::type);
        h= std::vector<int>(static_cast<unsigned long>(D.cols));
        obs=std::vector<obst>(static_cast<unsigned long>(D.cols)+1);
    }
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

void dt::Marie19::process1d(const cv::Mat &in)
{
    // Cf article (Marie19, Hirata96 ou Feizenszwalb12 pour l'algo)

    int rows=in.rows;
    int cols=in.cols;
    int maxD=rows+cols;

    // Légère modification pour parcourir l'image dans le bon sens
    // Gain non négligeable en temps de calcul

    auto m = in.ptr<uchar>(0);
    auto g = g_1d.ptr<int>(0);
    // Pour chaque première case de colonne
    for(int j=0;j<cols;j++){
        if(m[j]==0)
            g[j]=0;
        else
            g[j]=maxD;
    }

    // Même si on regarde les colonnes, on parcourt les lignes,
    // C'est beaucoup plus efficace
    for(int i=1;i<rows;i++){
        m = in.ptr<uchar>(i);
        g = g_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(m[j]==0)
                g[j]=0;
            else
                g[j]=g_1d.at<int>(i-1,j)+1;
        }
    }

    // Encore une fois, on regarde les colonnes, mais on parcourt les lignes
    for(int i=rows-2;i>=0;i--){
        g = g_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(g_1d.at<int>(i+1,j)<g[j]){
                g[j]=g_1d.at<int>(i+1,j)+1;
            }
        }
    }
}

void dt::Marie19::process1d_withP(const cv::Mat &in)
{
    // Cf article (Marie19, Hirata96 ou Feizenszwalb12 pour l'algo)

    int rows=in.rows;
    int cols=in.cols;
    int maxD=rows+cols;

    // Légère modification pour parcourir l'image dans le bon sens
    // Gain non négligeable en temps de calcul

    auto m = in.ptr<uchar>(0);
    auto g = g_1d.ptr<int>(0);
    auto p = p_1d.ptr<int>(rows-1);

    // Pour chaque première case de colonne
    for(int j=0;j<cols;j++){
        if(m[j]==0)
            g[j]=0;
        else
            g[j]=maxD;
    }

    // Même si on regarde les colonnes, on parcourt les lignes,
    // C'est beaucoup plus efficace
    for(int i=1;i<rows;i++){
        m = in.ptr<uchar>(i);
        g = g_1d.ptr<int>(i);
        for(int j=0;j<cols;j++){
            if(m[j]==0)
                g[j]=0;
            else
                g[j]=g_1d.at<int>(i-1,j)+1;
        }
    }

    // Calcul des projections 1d
    for(int j=0;j<cols;j++){
        if(g[j]>=maxD)
            p[j]=-1;
        else
            p[j]=rows-1-g[j];
    }

    // On étudie les colonnes, mas on parcourt les lignes
    // par soucis d'efficacité
    for(int i=rows-2;i>=0;i--){
        g = g_1d.ptr<int>(i);
        p = p_1d.ptr<int>(i);
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

void dt::Marie19::process2d(bool squared, cv::Mat &D)
{
    // Cf article (Marie19) pour l'algo

    // Variables de travail
    int rows=D.rows, cols=D.cols, maxD2=(rows+cols)*(rows+cols);
    int cur, w;

    // Bornes infinies pour convenance
    int infty = std::numeric_limits<int>::max();
    int minfty = std::numeric_limits<int>::min();


    // Pour chaque ligne
    for(int i=0; i < rows; i++) {

        // On positionne une parabole sur l'ensemble de la ligne
        auto q=obs.begin();
        q->s=0; q->t=minfty; (q+1)->t=infty; q->k=0; w=0;
        auto pG=g_1d.ptr<int>(i);

        // Optimisation 1 : On stocke la ligne dans un vecteur de travail
        for(unsigned long j=0;j<static_cast<unsigned long>(cols);j++)
            h[j]= pG[j]*pG[j];

        // Pour chaque case de la ligne courante
        for(unsigned long j=1;j<static_cast<unsigned long>(cols);j++) {
            // Optimisation 2 : si la cellule à une distance au bord infinie, on ne la considère pas
            cur=h[j];
            if(cur>=maxD2)
                continue;

            // Si
            // - D courant == D précédent
            // - Début de l'obstacle courant (w) < cellule courante
            if(cur==h[j-1]&&w<static_cast<int>(j)){
                // On passe à la cellule suivante
                j++;
                // Si cette cellule existe et qu'elle a également la même distance au bord
                if(static_cast<int>(j)<cols-1&&h[j]==cur){
                    // On injecte un segment
                    q++; q->t=static_cast<int>(j)-1; q->s=pG[j]; q->k=1;
                    // Et on fast forward jusqu'à la fin de celui-ci
                    while(static_cast<int>(j)<cols&&h[j]==cur)
                        j++;
                }
                // On ajoute la parabole de fin
                q++; q->s=static_cast<int>(j)-1;   q->k=0;  q->t=static_cast<int>(j)-1;   (q+1)->t=infty;

                // Si on est arrivé au bout de la ligne, on sort de la boucle
                if(static_cast<int>(j)>=cols)
                    break;
            }

            // Qu'on soit passé dans le if ou non, on cherche maintenant où insérer la parabole courante

            // On calcule l'intersection entre la parabole centrée en j et la parabole courante q (c'est forcément une parabole ici)
            w=static_cast<int>(1+SEP(q->s,static_cast<int>(j)));

            // Si cette intersection est hors zone, on passe à la cellule suivante
            if(w>=static_cast<int>(cols))
                continue;

            // Tant que cette intersection est inférieure au début de la zone d'influence de l'obstacle q
            while(w<=q->t){
                // On supprime l'obstacle courant
                q--;

                // Et on calcule l'intersection avec le nouvel obstacle courant
                if(!q->k) // Si c'est une parabole
                    w=static_cast<int>(1+SEP(q->s,static_cast<int>(j)));
                else // Si c'est un segment
                    w=static_cast<int>(1+SEP2(q->s,static_cast<int>(j)));
            }
            // On finit par ajouter la parabole qu'on vient de définir
            q++;   q->s=static_cast<int>(j);  q->t=w;   (q+1)->t=infty;   q->k=0;
        }

        // Reste à parcourir les obstacles pour assigner à chaque cellule sa distance au bord (intermédiaire)
        auto pD = D.ptr<double>(i);
        q=obs.begin();
        // Pour chaque case
        for(int j=0;j<cols;j++){
            // On se place sur le bon intervalle
            while((q+1)->t<=j)    q++;
            // On calcule D
            if(!q->k){ // a - si l'obstacle est une parabole
                if(squared)
                    pD[j]=F(q->s,j);
                else
                    pD[j]=sqrt(F(q->s,j));
            }
            else{      // b - si l'obstacle est un segment
                if(squared)
                    pD[j]=q->s*q->s;
                else
                    pD[j]=q->s;
            }
        }
    }
}

void dt::Marie19::process2d_withP(bool squared, cv::Mat &D, cv::Mat &P)
{
    // Cf article (Marie19) pour l'algo

    // Variables de travail
    int rows=D.rows, cols=D.cols, maxD2=(rows+cols)*(rows+cols);
    int cur, w;

    // Bornes infinies pour convenance
    int infty = std::numeric_limits<int>::max();
    int minfty = std::numeric_limits<int>::min();

    // Pour chaque ligne
    for(int i=0; i < rows; i++) {

        // On positionne une parabole sur l'ensemble de la ligne
        auto q=obs.begin();
        q->s=0; q->t=minfty; (q+1)->t=infty; q->k=0; w=0;
        auto pG=g_1d.ptr<int>(i);

        // Optimisation 1 : On stocke la ligne dans un vecteur de travail
        for(unsigned long j=0;j<static_cast<unsigned long>(cols);j++)
            h[j]= pG[j]*pG[j];

        // Pour chaque case de la ligne courante
        for(unsigned long j=1;j<static_cast<unsigned long>(cols);j++) {
            // Optimisation 2 : si la cellule à une distance au bord infinie, on ne la considère pas
            cur=h[j];
            if(cur>=maxD2)
                continue;

            // Si
            // - D courant == D précédent
            // - Début de l'obstacle courant (w) < cellule courante
            if(cur==h[j-1]&&w<static_cast<int>(j)){
                // On passe à la cellule suivante
                j++;
                // Si cette cellule existe et qu'elle a également la même distance au bord
                if(static_cast<int>(j)<cols-1&&h[j]==cur){
                    // On injecte un segment
                    q++; q->t=static_cast<int>(j)-1; q->s=pG[j]; q->k=1;
                    // Et on fast forward jusqu'à la fin de celui-ci
                    while(static_cast<int>(j)<cols&&h[j]==cur)
                        j++;
                }
                // On ajoute la parabole de fin
                q++; q->s=static_cast<int>(j)-1;   q->k=0;  q->t=static_cast<int>(j)-1;   (q+1)->t=infty;

                // Si on est arrivé au bout de la ligne, on sort de la boucle
                if(static_cast<int>(j)>=cols)
                    break;
            }

            // Qu'on soit passé dans le if ou non, on cherche maintenant où insérer la parabole courante

            // On calcule l'intersection entre la parabole centrée en j et la parabole courante q (c'est forcément une parabole ici)
            w=static_cast<int>(1+SEP(q->s,static_cast<int>(j)));

            // Si cette intersection est hors zone, on passe à la cellule suivante
            if(w>=static_cast<int>(cols))
                continue;

            // Tant que cette intersection est inférieure au début de la zone d'influence de l'obstacle q
            while(w<=q->t){
                // On supprime l'obstacle courant
                q--;

                // Et on calcule l'intersection avec le nouvel obstacle courant
                if(!q->k) // Si c'est une parabole
                    w=static_cast<int>(1+SEP(q->s,static_cast<int>(j)));
                else // Si c'est un segment
                    w=static_cast<int>(1+SEP2(q->s,static_cast<int>(j)));
            }
            // On finit par ajouter la parabole qu'on vient de définir
            q++;   q->s=static_cast<int>(j);  q->t=w;   (q+1)->t=infty;   q->k=0;
        }

        // Reste à parcourir les obstacles pour assigner à chaque cellule sa distance au bord (intermédiaire)
        auto pD = D.ptr<double>(i);
        auto pP = P.ptr<cv::Point>(i);
        auto pP1 = p_1d.ptr<int>(i);
        q=obs.begin();
        // Pour chaque case
        for(int j=0;j<cols;j++){
            // On se place sur le bon intervalle
            while((q+1)->t<=j)    q++;
            // On calcule D
            if(!q->k){ // a - si l'obstacle est une parabole
                if(squared)
                    pD[j]=F(q->s,j);
                else
                    pD[j]=sqrt(F(q->s,j));
                pP[j]=cv::Point(pP1[q->s],q->s);
            }
            else{      // b - si l'obstacle est un segment
                if(squared)
                    pD[j]=q->s*q->s;
                else
                    pD[j]=q->s;
                pP[j]=cv::Point(pP1[j],j);
            }
        }
    }
}

dt::Marie19_nd::Marie19_nd():AlgoNd("Marie 2019")
{
}

void dt::Marie19_nd::compute(const Matrice_ND<uchar> &in, Matrice_ND<double> &D, bool squared)
{
    // Si la matrice est vide...
    if(in.empty())
        return;

    // Si la matrice résultat est vide, on alloue la mémoire
    if(D.empty())
        D = Matrice_ND<double>(in.dims());

    // Bornes infinies pour convenance
    int infty = std::numeric_limits<int>::max();
    int minfty = std::numeric_limits<int>::min();

    int maxD = 0;
    for(int i:in.dims()){
        maxD+=i;
    }
    int maxD2 = maxD*maxD; // On borne la distance entre 2 cellules
    int w=0; // Indice de l'intersection entre 2 obstacles le long d'une ligne

    // On récupère les itérateurs de ligne suivant la direction 0
    std::vector<std::vector<unsigned long> > lignes;
    D.getLines(0,lignes);

    // Pour chaque ligne suivant la dimension 0
    for(std::vector<unsigned long> l:lignes){

        if(in.at(l[0])!=0) // S'il ne s'agit pas d'un point obstacle
            D[l[0]]=maxD;  // On initialise la distance de la première case à l'infini
        else               // Sinon
            D[l[0]]=0;     // On l'initialise à 0

        // Pour chaque case suivante
        for(unsigned long j=1;j<l.size();j++){
            if(in.at(l[j])!=0) // Si ce n'est pas un point obstacle
                D[l[j]]=1+D[l[j-1]]; // Sa distance au bord est 1 + celle de la cellule précédente
            else // Sinon
                D[l[j]]=0; // Sa distance au bord est égale à 0
        }

        // On parcourt maintenant la ligne dans l'autre sens
        for(int j=static_cast<int>(l.size())-2;j>=0;j--){
            // Si la distance au bord de la cellule à droite est inférieure à celle de la cellule courante
            if(D[l[static_cast<unsigned long>(j)+1]]<D[l[static_cast<unsigned long>(j)]])
                // Celle de la cellule courante devient égale à 1 + celle de la cellule à droite
                D[l[static_cast<unsigned long>(j)]]=1+D[l[static_cast<unsigned long>(j)+1]];
        }
        // Pour convenance, on met la carte des distances au carré (ne change rien au résultat final
        for(unsigned long j=0;j<l.size();j++){
            D[l[j]]*=D[l[j]];
        }
    }

    // Pour chacune des dimensions suivantes
    for(int i=1;i<in.nDim();i++){

        // On initialise les variables de travail
        unsigned long taille=static_cast<unsigned long>(D.dim(static_cast<unsigned long>(i)));
        std::vector<obst> obs(1+taille);
        std::vector<int> h(taille);
        int cur;

        // On récupère les itérateurs de ligne suivante la dimension courante
        D.getLines(i,lignes);

        // Pour chaque ligne de la dimension courante
        for(std::vector<unsigned long> l:lignes){

            // On positionne une parabole sur l'ensemble de la ligne
            auto q=obs.begin();
            q->s=0; q->t=minfty; (q+1)->t=infty; q->k=0;

            // Optimisation 1 : On stocke la ligne dans un vecteur de travail
            for(unsigned long j=0;j<taille;j++)
                h[j]= static_cast<int>(D[l[j]]);

            // Pour chaque cellule de la ligne courante
            for(unsigned long j=1;j<taille;j++) {
                // Optimisation 2 : si la cellule à une distance au bord infinie, on ne la considère pas
                cur=h[j];
                if(cur>=maxD2)
                    continue;

                // Si
                // - D courant == D précédent
                // - Début de l'obstacle courant (w) < cellule courante
                if(cur==h[j-1]&&w<static_cast<int>(j)){
                    // On passe à la cellule suivante
                    j++;
                    // Si cette cellule existe et qu'elle a également la même distance au bord
                    if(j<taille-1&&h[j]==cur){
                        // On injecte un segment
                        q++; q->t=static_cast<int>(j)-1; q->s=cur; q->k=1;
                        // Et on fast forward jusqu'à la fin de celui-ci
                        while(j<taille&&h[j]==cur)
                            j++;
                    }
                    // On ajoute la parabole de fin
                    q++; q->s=static_cast<int>(j)-1;   q->k=0;  q->t=static_cast<int>(j)-1;   (q+1)->t=infty;

                    // Si on est arrivé au bout de la ligne, on sort de la boucle
                    if(j>=taille)
                        break;

                }

                // Qu'on soit passé dans le if ou non, on cherche maintenant où insérer la parabole courante

                // On calcule l'intersection entre la parabole centrée en j et la parabole courante q (c'est forcément une parabole ici)
                w=1+SEP(q->s,static_cast<int>(j));

                // Si cette intersection est hors zone, on passe à la cellule suivante
                if(w>=static_cast<int>(taille))
                    continue;

                // Tant que cette intersection est inférieure au début de la zone d'influence de l'obstacle q
                while(w<=q->t){
                    // On supprime l'obstacle courant
                    q--;

                    // Et on calcule l'intersection avec le nouvel obstacle courant
                    if(!q->k) // Si c'est une parabole
                        w=1+SEP(q->s,static_cast<int>(j));
                    else // Si c'est un segment
                        w=1+SEP3(q->s,j);
                }
                // On finit par ajouter la parabole qu'on vient de définir
                q++;   q->s=static_cast<int>(j);  q->t=w;   (q+1)->t=infty;   q->k=0;
            }

            // Reste à parcourir les obstacles pour assigner à chaque cellule sa distance au bord (intermédiaire)
            q=obs.begin();
            // Pour chaque cellule
            for(unsigned long j=0;j<taille;j++){
                // On se place sur le bon intervalle
                while((q+1)->t<=static_cast<int>(j))    q++;
                // On calcule D
                if(!q->k){ // a - si l'obstacle est une parabole
                    D[l[j]]=F(q->s,static_cast<int>(j));
                }
                else{      // b - si l'obstacle est un segment
                    D[l[j]]=q->s;
                }
            }
        }
    }

    // Si on souhaite la distance euclidienne et non son carré
    if(!squared){
        // On applique sqrt à chaque cellule
        for(auto it(D.begin());it!=D.end();it++)
            *it=sqrt(*it);
    }
}

void dt::Marie19_nd::compute(const Matrice_ND<uchar> &in, Matrice_ND<double> &D, Matrice_ND<int > &P, bool squared)
{
    // Si la matrice est vide...
    if(in.empty())
        return;

    // Si la matrice résultat est vide, on alloue la mémoire
    if(D.empty())
        D = Matrice_ND<double>(in.dims());
    if(P.empty())
        P = Matrice_ND<int>(in.dims());

    // Bornes infinies pour convenance
    int infty = std::numeric_limits<int>::max();
    int minfty = std::numeric_limits<int>::min();

    int maxD = 0;
    for(int i:in.dims()){
        maxD+=i;
    }
    int maxD2 = maxD*maxD; // On borne la distance entre 2 cellules
    int w=0; // Indice de l'intersection entre 2 obstacles le long d'une ligne

    // On récupère les itérateurs de ligne suivant la direction 0
    std::vector<std::vector<unsigned long> > lignes;
    D.getLines(0,lignes);

    // Pour chaque ligne suivant la dimension 0
    for(std::vector<unsigned long> l:lignes){

        if(in.at(l[0])!=0) // S'il ne s'agit pas d'un point obstacle
            D[l[0]]=maxD;  // On initialise la distance de la première case à l'infini
        else{             // Sinon
            D[l[0]]=0;     // On l'initialise à 0
            P[l[0]]=static_cast<int>(l[0]);
        }

        // Pour chaque case suivante
        for(unsigned long j=1;j<l.size();j++){
            if(in.at(l[j])!=0){ // Si ce n'est pas un point obstacle
                D[l[j]]=1+D[l[j-1]]; // Sa distance au bord est 1 + celle de la cellule précédente
                P[l[j]]=P[l[j-1]];
            }
            else{ // Sinon
                D[l[j]]=0; // Sa distance au bord est égale à 0
                P[l[j]]=static_cast<int>(l[j]);
            }
        }

        // On parcourt maintenant la ligne dans l'autre sens
        for(int j=static_cast<int>(l.size())-2;j>=0;j--){
            // Si la distance au bord de la cellule à droite est inférieure à celle de la cellule courante
            if(D[l[static_cast<unsigned long>(j)+1]]<D[l[static_cast<unsigned long>(j)]]){
                // Celle de la cellule courante devient égale à 1 + celle de la cellule à droite
                D[l[static_cast<unsigned long>(j)]]=1+D[l[static_cast<unsigned long>(j)+1]];
                P[l[static_cast<unsigned long>(j)]]=P[l[static_cast<unsigned long>(j)+1]];
            }
        }
        // Pour convenance, on met la carte des distances au carré (ne change rien au résultat final
        for(unsigned long j=0;j<l.size();j++){
            if(in.at(l[j])!=0)
                D[l[j]]*=D[l[j]];
        }
    }

    // Pour chacune des dimensions suivantes
    for(int i=1;i<in.nDim();i++){

        // On initialise les variables de travail
        unsigned long taille=static_cast<unsigned long>(D.dim(static_cast<unsigned long>(i)));
        std::vector<obst> obs(1+taille);
        std::vector<int> h(taille);
        std::vector<int> p(taille);
        int cur;

        // On récupère les itérateurs de ligne suivante la dimension courante
        D.getLines(i,lignes);

        // Pour chaque ligne de la dimension courante
        int n=-1;
        for(std::vector<unsigned long> l:lignes){
            n++;



            // On positionne une parabole sur l'ensemble de la ligne
            auto q=obs.begin();
            q->s=0; q->t=minfty; (q+1)->t=infty; q->k=0; w=0;

            // Optimisation 1 : On stocke la ligne dans un vecteur de travail
            for(unsigned long j=0;j<taille;j++){
                h[j]= static_cast<int>(D[l[j]]);
                p[j]= P[l[j]];
            }

            // Pour chaque cellule de la ligne courante
            for(unsigned long j=1;j<taille;j++) {
                // Optimisation 2 : si la cellule à une distance au bord infinie, on ne la considère pas
                cur=h[j];
                if(cur>=maxD2)
                    continue;

                // Si
                // - D courant == D précédent
                // - Début de l'obstacle courant (w) < cellule courante
                if(cur==h[j-1]&&w<static_cast<int>(j)){
                    // On passe à la cellule suivante
                    j++;
                    // Si cette cellule existe et qu'elle a également la même distance au bord
                    if(j<taille-1&&h[j]==cur){
                        // On injecte un segment
                        q++; q->t=static_cast<int>(j)-1; q->s=cur; q->k=1;
                        // Et on fast forward jusqu'à la fin de celui-ci
                        while(j<taille&&h[j]==cur)
                            j++;
                    }
                    // On ajoute la parabole de fin
                    q++; q->s=static_cast<int>(j)-1;   q->k=0;  q->t=static_cast<int>(j)-1;   (q+1)->t=infty;

                    // Si on est arrivé au bout de la ligne, on sort de la boucle
                    if(j>=taille)
                        break;

                }

                // Qu'on soit passé dans le if ou non, on cherche maintenant où insérer la parabole courante

                // On calcule l'intersection entre la parabole centrée en j et la parabole courante q (c'est forcément une parabole ici)
                w=static_cast<int>(1+SEP(q->s,static_cast<int>(j)));

                // Si cette intersection est hors zone, on passe à la cellule suivante
                if(w>=static_cast<int>(taille))
                    continue;

                // Tant que cette intersection est inférieure au début de la zone d'influence de l'obstacle q
                while(w<=q->t){
                    // On supprime l'obstacle courant
                    q--;

                    // Et on calcule l'intersection avec le nouvel obstacle courant
                    if(!q->k) // Si c'est une parabole
                        w=static_cast<int>(1+SEP(q->s,static_cast<int>(j)));
                    else // Si c'est un segment
                        w=static_cast<int>(1+SEP3(q->s,j));
                }
                // On finit par ajouter la parabole qu'on vient de définir
                q++;   q->s=static_cast<int>(j);  q->t=w;   q->k=0;   (q+1)->t=infty;
            }

            // Reste à parcourir les obstacles pour assigner à chaque cellule sa distance au bord (intermédiaire)
            q=obs.begin();
            // Pour chaque cellule
            for(unsigned long j=0;j<taille;j++){
                // On se place sur le bon intervalle
                while((q+1)->t<=static_cast<int>(j))    q++;
                // On calcule D
                if(!q->k){ // a - si l'obstacle est une parabole
                    D[l[j]]=F(q->s,static_cast<int>(j));
                    P[l[j]]=p[static_cast<unsigned long>(q->s)];
                }
                else{      // b - si l'obstacle est un segment
                    D[l[j]]=q->s;
                }
            }
        }
    }

    // Si on souhaite la distance euclidienne et non son carré
    if(!squared){
        // On applique sqrt à chaque cellule
        for(auto it(D.begin());it!=D.end();it++)
            *it=sqrt(*it);
    }
}
