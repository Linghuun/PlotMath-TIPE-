#ifndef QUAD
#define QUAD
#include "point.h"

class Polynome;
class Camera;
class Quadtree;

typedef struct int_pile{int val; int_pile* next;} int_pile;
typedef struct quad_courb { Quadtree* quad;
                            double courb;} quad_courb;
quad_courb* quad_courb_creer(Quadtree* q, double c);

class Quadtree{
    public:
        Quadtree(int subdiv, int nb_max_point, Polynome* poly);
        ~Quadtree();
        void copier(Quadtree* quadtree);

        bool appartient(std::shared_ptr<Point> point);
        void add_point(std::shared_ptr<Point> point); 
        void push_vector(std::vector<std::shared_ptr<Point>>* points); 
        void divide_space();
        int get_point_count();

        std::shared_ptr<Point> get_center(); 
        bool is_divided();
        void nettoyer_double(double x_min, double x_max, double y_min, double y_max, int n, double precision);
        void lier_points(std::vector<std::shared_ptr<Point>>* non_lies, std::vector<int>* partition, int& partition_id);
        int zone_non_vide(double x_min, double x_max, double y_min, double y_max);
        bool zone_recherche(double& x_min, double& x_max, double& y_min, double& y_max);
        double courbure_moy_quadtree();
        
        void augmenter_qualite_visible(Camera* camera, int nb_iter, double precision, int version, Quadtree* self,
                                        bool qualite, std::vector<quad_courb*>* tableau_quad_courb);
        void augmenter_qualite_V1_ajout(int nb_iter, double precision);
        void augmenter_qualite_V2_creation(int nb_iter, double precision, double facteur);

        Quadtree** quadtrees;
        std::vector<std::shared_ptr<Point>>* points;
        Polynome* forme;
        double min_x, min_y, max_x, max_y;
        bool divided;
        std::shared_ptr<Point> center; 
        int subdivision;
        int nb_pts_par_quadtree;
        int point_counter;
};
#endif
