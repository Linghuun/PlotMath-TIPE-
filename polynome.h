#ifndef POLY
#define POLY

#include <vector>
#include "point.h"
class Quadtree;

class Terme{
    public:
        Terme(int pow_x, int pow_y, double coeff);

        Terme* deriver_x();
        Terme* deriver_y();
        double get_coeff();
        int pow_x;
        int pow_y;
        double coeff;
};

class Polynome{
    public:
        Polynome();
        ~Polynome();

        std::vector<Terme*>* termes;

        void ajouter_terme(Terme* terme);

        Polynome* deriver_x();
        Polynome* deriver_y();
        Polynome* evaluer_x(double x);
        Polynome* evaluer_y(double y);
        double evaluer_xy(double x, double y);
        double evaluer_y_valeur(double y);
        double evaluer_x_valeur(double x);
        double valeur();
        void print(char* name);

        Polynome* multiplier_polynome(Polynome* p_2);
        Polynome* additionner_polynome(Polynome* p_2);
        void calculer_courbure(Quadtree* quadtree);
        double calculer_courbure_xy(double x, double y);
        
        std::vector<Point*>* trouver_racine_y(double x, double y_min, double y_max, int nb_iter, int nb_pts, double precision, double multiplicateur);
        std::vector<Point*>* trouver_racine_x(double y, double x_min, double x_max, int nb_iter, int nb_pts, double precision, double multiplicateur);
        double trouver_racine_newton_y(int nb_iter, double y_0, Polynome* poly_deriv, double precision);
        double trouver_racine_newton_x(int nb_iter, double x_0, Polynome* poly_deriv, double precision);

        Quadtree* generer_quadtree_poly(double x_max, double x_min, double y_max, double y_min, int nb_iter, int nb_pts,
                                                            int nb_pts_par_quadtree, double precision, double multiplicateur);

    protected:
        
};

#endif