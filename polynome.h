#ifndef POLY
#define POLY

#include <vector>
#include <memory>
#include "point.h"
class Quadtree;

class Terme {
    public:
        Terme(int pow_x, int pow_y, double coeff);
        std::unique_ptr<Terme> deriver_x_terme();
        std::unique_ptr<Terme> deriver_y_terme();
        
        int pow_x;
        int pow_y;
        double coeff;
};

class Polynome {
    public:
        Polynome();
        ~Polynome();
        
        std::vector<std::unique_ptr<Terme>>* termes;
        void ajouter_terme(std::unique_ptr<Terme> terme);

        Polynome* deriver_x();
        Polynome* deriver_y();
        Polynome* evaluer_x(double x);
        Polynome* evaluer_y(double y);
        double evaluer_xy(double x, double y);
        double evaluer_y_valeur(double y);
        double evaluer_x_valeur(double x);

        void print(char* name);
        
        std::vector<std::shared_ptr<Point>>* trouver_racine_y(double x, double y_min, double y_max, int nb_iter, int nb_pts, double precision);
        std::vector<std::shared_ptr<Point>>* trouver_racine_x(double y, double x_min, double x_max, int nb_iter, int nb_pts, double precision);
        double trouver_racine_newton_y(int nb_iter, double y_0, Polynome* poly_deriv, double precision);
        double trouver_racine_newton_x(int nb_iter, double x_0, Polynome* poly_deriv, double precision);

        Quadtree* generer_quadtree_poly(double x_max, double x_min, double y_max, double y_min, int nb_iter, int nb_pts, int nb_pts_par_quadtree,
                                            int nb_section_doublons, double precision, Polynome* self, std::vector<std::shared_ptr<Point>>* non_lies);
};

#endif
