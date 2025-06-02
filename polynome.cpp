#include "polynome.h"
#include "quadtree.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <memory>
#include <cstdlib>

// Implementation Terme
Terme::Terme(int pow_x_, int pow_y_, double coeff_) {
    pow_x = pow_x_;
    pow_y = pow_y_;
    coeff = coeff_;
}

std::unique_ptr<Terme> Terme::deriver_x_terme() {
    if (pow_x == 0) {
        return std::unique_ptr<Terme>(new Terme(0, 0, 0));
    }
    return std::unique_ptr<Terme>(new Terme(pow_x - 1, pow_y, coeff * pow_x));
}

std::unique_ptr<Terme> Terme::deriver_y_terme() {
    if (pow_y == 0) {
        return std::unique_ptr<Terme>(new Terme(0, 0, 0));
    }
    return std::unique_ptr<Terme>(new Terme(pow_x, pow_y - 1, coeff * pow_y));
}
// Fin implementation Terme

// Implementation Polynome
Polynome::Polynome() {
    termes = new std::vector<std::unique_ptr<Terme>>;
}

Polynome::~Polynome() {
    delete termes;
}

void Polynome::print(char* name) {
    printf("%s = ", name);
    for (size_t i = 0; i < termes->size(); i++) {
        auto& terme = (*termes)[i];
        printf("%f X^%d Y^%d + ", terme->coeff, terme->pow_x, terme->pow_y);
    }
    printf("\n");
}

void Polynome::ajouter_terme(std::unique_ptr<Terme> terme) {
    if (terme->coeff != 0) {
        termes->push_back(std::move(terme));
    }
}

Polynome* Polynome::deriver_x() {
    Polynome* poly = new Polynome();
    for (size_t i = 0; i < termes->size(); i++) {
        auto& t = (*termes)[i];
        auto t_derive = t->deriver_x_terme();
        poly->ajouter_terme(std::move(t_derive));
    }
    return poly;
}

Polynome* Polynome::deriver_y() {
    Polynome* poly = new Polynome();
    for (size_t i = 0; i < termes->size(); i++) {
        auto& t = (*termes)[i];
        auto t_derive = t->deriver_y_terme();
        poly->ajouter_terme(std::move(t_derive));
    }
    return poly;
}

Polynome* Polynome::evaluer_x(double val_x) {
    Polynome* p = new Polynome();
    size_t nb_termes = termes->size();
    for (size_t i = 0; i < nb_termes; i++) {
        std::unique_ptr<Terme> terme(new Terme(
            0, 
            (*termes)[i]->pow_y, 
            (*termes)[i]->coeff * std::pow(val_x, (*termes)[i]->pow_x)
        ));
        p->ajouter_terme(std::move(terme));
    }
    return p;
}

Polynome* Polynome::evaluer_y(double val_y) {
    Polynome* p = new Polynome();
    size_t nb_termes = termes->size();
    for (size_t i = 0; i < nb_termes; i++) {
        std::unique_ptr<Terme> terme(new Terme(
            (*termes)[i]->pow_x, 
            0, 
            (*termes)[i]->coeff * std::pow(val_y, (*termes)[i]->pow_y)
        ));
        p->ajouter_terme(std::move(terme));
    }
    return p;
}

double Polynome::evaluer_xy(double x, double y) {
    double val = 0;
    size_t n = termes->size();
    for (size_t i = 0; i < n; i++) {
        val += (*termes)[i]->coeff * std::pow(y, (*termes)[i]->pow_y) * std::pow(x, (*termes)[i]->pow_x);
    }
    return val;
}

double Polynome::evaluer_y_valeur(double val_y) {
    double val = 0;
    size_t n = termes->size();
    for (size_t i = 0; i < n; i++) {
        if (!(*termes)[i]) {
            std::cerr << "Terme[" << i << "] est nul" << std::endl;
            continue; // ou return NAN selon le contexte
        }

        val += (*termes)[i]->coeff * std::pow(val_y, (*termes)[i]->pow_y);
    }
    return val;
}

double Polynome::evaluer_x_valeur(double val_x) {
    double val = 0;
    size_t n = termes->size();
    for (size_t i = 0; i < n; i++) {
        if (!(*termes)[i]) {
            std::cerr << "Terme[" << i << "] est nul" << std::endl;
            continue; // ou return NAN selon le contexte
        }
        val += (*termes)[i]->coeff * std::pow(val_x, (*termes)[i]->pow_x);
    }
    return val;
}

double Polynome::trouver_racine_newton_y(int nb_iter, double y_0, Polynome* poly_deriv, double precision) {
    double y_n = y_0;
    double y_n1 = y_0;
    for (int i = 0; i < nb_iter; i++) {
        y_n = y_n1;
        double f_yn = evaluer_y_valeur(y_n);
        double f_deriv_yn = poly_deriv->evaluer_y_valeur(y_n);
        if (f_deriv_yn == 0) {
            return NAN; // Evite la division par zéro
        }
        y_n1 = y_n - (f_yn / f_deriv_yn);

        if (std::abs(y_n1 - y_n) < precision) {
            return y_n1; // convergé, on peut sortir plus tôt
        }
    }
    return NAN;
}

double Polynome::trouver_racine_newton_x(int nb_iter, double x_0, Polynome* poly_deriv, double precision) {
    double x_n = x_0;
    double x_n1 = x_0;
    for (int i = 0; i < nb_iter; i++) {
        x_n = x_n1;
        double f_xn = evaluer_x_valeur(x_n);
        double f_deriv_xn = poly_deriv->evaluer_x_valeur(x_n);
        if (f_deriv_xn == 0) {
            return NAN; // Evite la division par zéro
        }
        x_n1 = x_n - (f_xn / f_deriv_xn);
    }

    if (std::abs(x_n - x_n1) < precision) {
        return x_n1;
    }
    return NAN;
}

std::vector<std::shared_ptr<Point>>* Polynome::trouver_racine_y(double x, double y_min, double y_max,
                                                                 int nb_iter, int nb_pts, double precision) {
    auto racines = new std::vector<std::shared_ptr<Point>>;
    std::unique_ptr<Polynome> p_eval(evaluer_x(x));
    std::unique_ptr<Polynome> p_deriv(p_eval->deriver_y());

    double pas = (y_max - y_min) / nb_pts;
    for (double y = y_min; y <= y_max; y += pas) {
        double racine = p_eval->trouver_racine_newton_y(nb_iter, y, p_deriv.get(), precision);
        if (racine > y_min && racine < y_max) {
            racines->push_back(std::make_shared<Point>(x, racine));
        }
    }
    return racines;
}

std::vector<std::shared_ptr<Point>>* Polynome::trouver_racine_x(double y, double x_min, double x_max,
                                                                 int nb_iter, int nb_pts, double precision) {
    auto racines = new std::vector<std::shared_ptr<Point>>;
    std::unique_ptr<Polynome> p_eval(evaluer_y(y));
    std::unique_ptr<Polynome> p_deriv(p_eval->deriver_x());

    double pas = (x_max - x_min) / nb_pts;
    for (double x = x_min; x <= x_max; x += pas) {
        double racine = p_eval->trouver_racine_newton_x(nb_iter, x, p_deriv.get(), precision);
        if (racine > x_min && racine < x_max) {
            racines->push_back(std::make_shared<Point>(racine, y));
        }
    }
    return racines;
}

Quadtree* Polynome::generer_quadtree_poly(double x_max, double x_min, double y_max, double y_min,
                                           int nb_iter, int nb_pts, int nb_pts_par_quadtree,
                                           int nb_section_doublons, double precision,
                                           Polynome* self, std::vector<std::shared_ptr<Point>>* non_lies){
    Quadtree* quadtree = new Quadtree(4, nb_pts_par_quadtree, self);
    quadtree->max_x = x_max;
    quadtree->min_x = x_min; 
    quadtree->max_y = y_max; 
    quadtree->min_y = y_min; 
    quadtree->center->set_x((x_min + x_max)/2);
    quadtree->center->set_y((y_min + y_max)/2);

    double pas_y = (y_max - y_min) / nb_pts;
    for (double y = y_min; y <= y_max; y += pas_y) {
        auto racines = trouver_racine_x(y, x_min, x_max, nb_iter, nb_pts, precision);
        quadtree->push_vector(racines);
        delete racines;
    }

    double pas_x = (x_max - x_min) / nb_pts;
    for (double x = x_min; x <= x_max; x += pas_x) {
        auto racines = trouver_racine_y(x, y_min, y_max, nb_iter, nb_pts, precision);
        quadtree->push_vector(racines);
        delete racines;
    }

    quadtree->nettoyer_double(x_min, x_max, y_min, y_max, nb_section_doublons, precision*0.1);

    quadtree->divide_space();
    
    auto partition = new std::vector<int>;
    int partition_id = 0;
    quadtree->lier_points(non_lies, partition, partition_id);
    delete partition;

    return quadtree;
}