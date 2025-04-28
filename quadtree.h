#ifndef QUAD
#define QUAD

#include <vector>
#include "point.h"

typedef struct int_pile{int val; int_pile* next;} int_pile;

class Polynome;

class Quadtree{
    public:
        Quadtree(int subdiv, int nb_max_point);
        ~Quadtree();

        void add_point(Point* point);
        void push_vector(std::vector<Point*>* points);
        void divide_space();
        int get_point_count();
        Point* get_point(int i);
        Point* get_center();
        bool is_divided();
        void nettoyer_double(double precision);
        void lier_points(std::vector<Point*>* non_lies, double pas_max);

        double calculer_courbure_xy(double x, double y);

        Quadtree** quadtrees;
        std::vector<Point*>* points;
        Polynome** courbure;
        bool courbure_init;

    protected:
        bool divided;
        Point* center;
        double min_x, min_y, max_x, max_y;
        int subdivision;
        int nb_pts_par_quadtree;
        int point_counter;
        
    
};

#endif