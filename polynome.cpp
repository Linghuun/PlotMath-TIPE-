#include "polynome.h"
#include "quadtree.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>

Terme::Terme(int pow_x_, int pow_y_, double coeff_){
    pow_x = pow_x_;
    pow_y = pow_y_;
    coeff = coeff_;
}

Terme* Terme::deriver_x(){
    if (pow_x == 0){
        return new Terme(0, 0, 0);
    }
    return new Terme(pow_x-1, pow_y, coeff*pow_x);
}

Terme* Terme::deriver_y(){
    if (pow_y == 0){
        return new Terme(0, 0, 0);
    }
    return new Terme(pow_x, pow_y-1, coeff*pow_y);
}

double Terme::get_coeff(){
    return coeff;
}


Polynome::Polynome(){
    termes = new std::vector<Terme*>;
}

Polynome::~Polynome(){
    for (size_t i = 0; i < termes->size(); i++){
        delete (*termes)[i];
    }
    delete termes;
}

void Polynome::print(char* name){
    printf("%s = ", name);
    for (size_t i = 0; i < termes->size(); i++){
        Terme* terme = (*termes)[i];
        printf("%f X^%d Y^%d + ", terme->coeff, terme->pow_x, terme->pow_y);
    }
    printf("\n");
}

void Polynome::ajouter_terme(Terme* terme){
    if (terme->get_coeff() != 0){
        termes->push_back(terme);
    }
    else{
        delete terme;
    }
}

Polynome* Polynome::deriver_x(){
    Polynome* poly = new Polynome();
    size_t nb_termes = termes->size();
    for (size_t i  = 0; i < nb_termes; i++){
        poly->ajouter_terme((*termes)[i]->deriver_x());
    }
    return poly;
}

Polynome* Polynome::deriver_y(){
    Polynome* poly = new Polynome();
    size_t nb_termes = termes->size();
    for (size_t i  = 0; i < nb_termes; i++){
        poly->ajouter_terme((*termes)[i]->deriver_y());
    }
    return poly;
}

double power(double x, int n){
    if(n == 0){
        return 1;
    }
    if(n%2 == 0){
        return power(x*x, n/2);
    }
    return x*power(x*x, n/2);
}

Polynome* Polynome::evaluer_x(double val_x){
    Polynome* p = new Polynome();
    size_t nb_termes = termes->size();
    for(size_t i = 0; i < nb_termes; i++){
        Terme* terme = new Terme(0, 0, 0);
        terme->pow_x = 0;
        terme->pow_y = (*termes)[i]->pow_y;
        terme->coeff = (*termes)[i]->coeff * power(val_x, (*termes)[i]->pow_x);
        p->ajouter_terme(terme);
    }
    return p;
}


Polynome* Polynome::evaluer_y(double val_y){
    Polynome* p = new Polynome();
    size_t nb_termes = termes->size();
    for(size_t i = 0; i < nb_termes; i++){ 
        Terme* terme = new Terme(0, 0, 0);
        terme->pow_y = 0;
        terme->pow_x = (*termes)[i]->pow_x;
        terme->coeff = (*termes)[i]->coeff * power(val_y, (*termes)[i]->pow_y);
        p->ajouter_terme(terme);
    }
    return p;
}

double Polynome::evaluer_xy(double x, double y){
    double val = 0;
    size_t n = termes->size();
    for(size_t i = 0; i < n; i++){
        val = val + ((*termes)[i]->coeff * power(y, (*termes)[i]->pow_y) * power(x, (*termes)[i]->pow_x) );
    }
    return val;
}

double Polynome::evaluer_y_valeur(double val_y){
    double val = 0;
    size_t n = termes->size();
    for (size_t i = 0; i < n; i++){
        val = val + ((*termes)[i]->coeff * power(val_y, (*termes)[i]->pow_y));
    }
    return val;
}

double Polynome::evaluer_x_valeur(double val_x){
    double val = 0;
    size_t n = termes->size();
    for (size_t i = 0; i < n; i++){
        val = val + ((*termes)[i]->coeff * power(val_x, (*termes)[i]->pow_x));
    }
    return val;
}

double absolute(double val){
    if (val < 0){
        return -val;
    }
    return val;
}
double Polynome::trouver_racine_newton_y(int nb_iter, double y_0, Polynome* poly_deriv, double precision){
    double y_n = y_0;
    double y_n1 = y_0;
    for (int i = 0; i < nb_iter; i++){
        y_n = y_n1;
        double f_yn = evaluer_y_valeur(y_n);
        double f_deriv_yn = poly_deriv->evaluer_y_valeur(y_n);
        y_n1 = y_n - (f_yn/f_deriv_yn);
    }

    if (absolute(y_n-y_n1) < precision){
        return y_n1; //On ne renvoie la racine que si elle a convergé
    }
    return NAN;
}

double Polynome::trouver_racine_newton_x(int nb_iter, double x_0, Polynome* poly_deriv, double precision){
    double x_n = x_0;
    double x_n1 = x_0;
    for (int i = 0; i < nb_iter; i++){
        x_n = x_n1;
        double f_xn = evaluer_x_valeur(x_n);
        double f_deriv_xn = poly_deriv->evaluer_x_valeur(x_n);
        x_n1 = x_n - (f_xn/f_deriv_xn);
    }

    if (absolute(x_n-x_n1) < precision){
        return x_n1; //On ne renvoie la racine que si elle a convergé
    }
    return NAN;
}


std::vector<Point*>* Polynome::trouver_racine_y(double x, double y_min, double y_max, int nb_iter, int nb_pts, double precision, double multiplicateur){
    std::vector<Point*>* racines = new std::vector<Point*>;
    Polynome* p_eval = evaluer_x(x);
    Polynome* p_deriv = p_eval->deriver_y();

    double pas = (y_max-y_min)/nb_pts;
    for (double y = y_min; y <= y_max; y = y+pas){
        double racine = p_eval->trouver_racine_newton_y(nb_iter, y, p_deriv, precision);
        if (racine < y_max && racine > y_min){ // NAN est ignoré ici.
            // bool nouvelle_racine = true;
            // for(size_t i = 0; i < racines->size(); i ++){
            //     if (absolute((*racines)[i]->y() - racine*multiplicateur) < precision*multiplicateur){
            //         nouvelle_racine = false;
            //     }
            // }
            // if (nouvelle_racine){
            Point* point = new Point(x*multiplicateur, racine*multiplicateur);
            racines->push_back(point);
            // }
        }
    }
    delete p_eval;
    delete p_deriv;
    return racines;
}

std::vector<Point*>* Polynome::trouver_racine_x(double y, double x_min, double x_max, int nb_iter, int nb_pts, double precision, double multiplicateur){
    std::vector<Point*>* racines = new std::vector<Point*>;
    Polynome* p_eval = evaluer_y(y);
    Polynome* p_deriv = p_eval->deriver_x();

    double pas = (x_max-x_min)/nb_pts;
    for (double x = x_min; x <= x_max; x = x+pas){
        double racine = p_eval->trouver_racine_newton_x(nb_iter, x, p_deriv, precision);
        if (racine < x_max && racine > x_min){ // NAN ne passera pas ce test
            // bool nouvelle_racine = true;
            // for(size_t i = 0; i < racines->size(); i ++){
            //     if (absolute((*racines)[i]->y() - racine*multiplicateur) < precision*multiplicateur){
            //         nouvelle_racine = false;
            //     }
            // }
            // if (nouvelle_racine){
            Point* point = new Point(racine*multiplicateur, y*multiplicateur);
            racines->push_back(point);

            // }
        }
    }
    delete p_eval;
    delete p_deriv;
    return racines;
}

Quadtree* Polynome::generer_quadtree_poly(double x_max, double x_min, double y_max, double y_min, int nb_iter, int nb_pts, int nb_pts_par_quadtree, double precision, double multiplicateur){
    Quadtree* quadtree = new Quadtree(4, nb_pts_par_quadtree);

    double pas_y = (y_max-y_min)/nb_pts;
    for (double y = y_min; y <= y_max; y = y + pas_y){
        std::vector<Point*>* racines = trouver_racine_x(y, x_min, x_max, nb_iter, nb_pts, precision, multiplicateur);
        quadtree->push_vector(racines);
        delete racines;
    }

    double pas_x = (x_max-x_min)/nb_pts;
    for (double x = x_min; x <= x_max; x = x + pas_x){
        std::vector<Point*>* racines = trouver_racine_y(x, y_min, y_max, nb_iter, nb_pts, precision, multiplicateur);
        quadtree->push_vector(racines);
        delete racines;
    }

    quadtree->nettoyer_double(precision*multiplicateur);
    quadtree->divide_space();

    std::vector<Point*>* non_lies = new std::vector<Point*>;
    double pas_max = 1.5*(pas_y+pas_x);

    quadtree->lier_points(non_lies, pas_max*multiplicateur);
    delete non_lies;
    calculer_courbure(quadtree);
    return quadtree;
}

Polynome* Polynome::multiplier_polynome(Polynome* p_2){
    Polynome* nouveau_poly = new Polynome();
    for(size_t i = 0; i < termes->size(); i++){
        for(size_t j = 0; j < p_2->termes->size(); j++){
            Terme* t1 = (*termes)[i];
            Terme* t2 = (*p_2->termes)[j];

            Terme* nouveau_terme = new Terme(0, 0, 0);
            nouveau_terme->pow_x = t1->pow_x + t2->pow_x;
            nouveau_terme->pow_y = t1->pow_y + t2->pow_y;
            nouveau_terme->coeff = t1->coeff * t2->coeff;

            nouveau_poly->ajouter_terme(nouveau_terme);
        }
    }
    
    return nouveau_poly;
}

Polynome* Polynome::additionner_polynome(Polynome* p_2){
    Polynome* nouveau_poly = new Polynome();
    for(size_t i = 0; i < termes->size(); i++){
        Terme* t = (*termes)[i];

        Terme* nouveau_terme = new Terme(0, 0, 0);
        nouveau_terme->pow_x = t->pow_x;
        nouveau_terme->pow_y = t->pow_y;
        nouveau_terme->coeff = t->coeff;

        nouveau_poly->ajouter_terme(nouveau_terme);
    }

    for(size_t i = 0; i < p_2->termes->size(); i++){
        Terme* t = (*p_2->termes)[i];
        
        Terme* nouveau_terme = new Terme(0, 0, 0);
        nouveau_terme->pow_x = t->pow_x;
        nouveau_terme->pow_y = t->pow_y;
        nouveau_terme->coeff = t->coeff;

        nouveau_poly->ajouter_terme(nouveau_terme);
    }
    return nouveau_poly;
}

void Polynome::calculer_courbure(Quadtree* quadtree){
    quadtree->courbure_init = true;
    quadtree->courbure = (Polynome**)malloc(sizeof(Polynome*)*5);

    quadtree->courbure[0] = deriver_y();
    quadtree->courbure[1] = deriver_x();
    quadtree->courbure[2] = quadtree->courbure[0]->deriver_y();
    quadtree->courbure[3] = quadtree->courbure[1]->deriver_x();
    quadtree->courbure[4] = quadtree->courbure[1]->deriver_y();   
}

