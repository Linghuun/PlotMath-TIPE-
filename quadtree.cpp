#include "quadtree.h"
#include "polynome.h"
#include "camera.h"

#include <cstdlib>
#include <iostream>
#include <float.h>
#include <algorithm>
#include <cmath>
#include <vector>

Quadtree::Quadtree(int subdiv, int nb_pts_max, Polynome* poly){
    quadtrees = (Quadtree**)malloc(sizeof(Quadtree*)*subdiv);

    divided = false;
    max_x = -DBL_MAX;
    max_y = -DBL_MAX;
    min_x = DBL_MAX;
    min_y = DBL_MAX;
    center = std::make_shared<Point>(0, 0);
    point_counter = 0;
    subdivision = subdiv;
    nb_pts_par_quadtree = nb_pts_max;
    points = new std::vector<std::shared_ptr<Point>>();
    forme = poly;
}

Quadtree::~Quadtree(){
    if (divided){
        for (int i = 0; i < subdivision; i++){
            delete quadtrees[i];
        }
    }
    else{
        if (points != nullptr){
            points->clear();
        }
    }
    free(quadtrees);
    delete points;
}

void Quadtree::add_point(std::shared_ptr<Point> point){
    point_counter++;
    points->push_back(point);
}

void Quadtree::push_vector(std::vector<std::shared_ptr<Point>>* points){
    for (size_t i = 0; i < points->size(); i++){
        add_point((*points)[i]);
    }
}

void Quadtree::divide_space(){
    if (point_counter >= nb_pts_par_quadtree){
        divided = true;
        for (int i = 0; i < subdivision; i++){
            quadtrees[i] = new Quadtree(subdivision, nb_pts_par_quadtree, forme);
        }
        quadtrees[0]->min_x = min_x;
        quadtrees[0]->max_y = max_y;
        quadtrees[0]->max_x = center->x();
        quadtrees[0]->min_y = center->y();
        quadtrees[0]->center->set_x((quadtrees[0]->min_x + quadtrees[0]->max_x)/2);
        quadtrees[0]->center->set_y((quadtrees[0]->min_y + quadtrees[0]->max_y)/2);


        quadtrees[1]->min_x = center->x();
        quadtrees[1]->max_y = max_y;
        quadtrees[1]->max_x = max_x;
        quadtrees[1]->min_y = center->y();
        quadtrees[1]->center->set_x((quadtrees[1]->min_x + quadtrees[1]->max_x)/2);
        quadtrees[1]->center->set_y((quadtrees[1]->min_y + quadtrees[1]->max_y)/2);

        quadtrees[2]->min_x = min_x;
        quadtrees[2]->max_y = center->y();
        quadtrees[2]->max_x = center->x();
        quadtrees[2]->min_y = min_y;
        quadtrees[2]->center->set_x((quadtrees[2]->min_x + quadtrees[2]->max_x)/2);
        quadtrees[2]->center->set_y((quadtrees[2]->min_y + quadtrees[2]->max_y)/2);

        quadtrees[3]->min_x = center->x();
        quadtrees[3]->max_y = center->y();
        quadtrees[3]->max_x = max_x;
        quadtrees[3]->min_y = min_y;
        quadtrees[3]->center->set_x((quadtrees[3]->min_x + quadtrees[3]->max_x)/2);
        quadtrees[3]->center->set_y((quadtrees[3]->min_y + quadtrees[3]->max_y)/2);

        for (size_t i = 0; i < points->size(); i++){
            if ((*points)[i]->x() < center->x()){
                if ((*points)[i]->y() < center->y()){
                    quadtrees[2]->add_point((*points)[i]);
                }
                else{
                    quadtrees[0]->add_point((*points)[i]);
                }
            }
            else{
                if ((*points)[i]->y() > center->y()){
                    quadtrees[1]->add_point((*points)[i]);
                }
                else{
                    quadtrees[3]->add_point((*points)[i]);
                }
            }
        }
        points->clear();
        quadtrees[0]->divide_space();
        quadtrees[1]->divide_space();
        quadtrees[2]->divide_space();
        quadtrees[3]->divide_space();
    }
}

bool doublon(int i, int j, std::vector<std::shared_ptr<Point>>* points, double precision){
    if (std::abs((*points)[i]->x() - (*points)[j]->x()) < precision){
        if (std::abs((*points)[i]->y() - (*points)[j]->y()) < precision){
            return true;
        }
    }
    return false;
}

void Quadtree::nettoyer_double(double x_min, double x_max, double y_min, double y_max, int nb_section_doublons, double precision){
    int n = nb_section_doublons;
    double section_x = std::fmax(std::abs(x_min), std::abs(x_max))/n;
    double section_y = std::fmax(std::abs(y_min), std::abs(y_max))/n;

    // Création de la grille
    std::vector<int>*** grille = (std::vector<int>***)malloc(sizeof(std::vector<int>**)*n);
    for (int i = 0; i < n; i++){
        grille[i] = (std::vector<int>**)malloc(sizeof(std::vector<int>*) * n);
        for (int j = 0; j < n; j++){
            grille[i][j] = new std::vector<int>;
        }
    }

    bool* a_supprimer = (bool*)malloc(sizeof(bool)*point_counter);
    for (int i = 0; i < point_counter; i++){
        a_supprimer[i] = false;
    }

    // Remplissage
    for (int i = 0; i < point_counter; i++){
        std::shared_ptr<Point> point = (*points)[i];

        //On peut avoir x = nb_sections_doublons si |x| = max(|x_min|,|x_max|)
        int x = std::min((int)(std::abs(point->x()) / section_x), nb_section_doublons-1);
        int y = std::min((int)(std::abs(point->y()) / section_y), nb_section_doublons-1);
        grille[x][y]->push_back(i);
    }

    // Comparaison
    for (int x = 0; x < n; x++){
        for (int y = 0; y < n; y++){
            if (grille[x][y] != NULL){
                std::vector<int>* pile = grille[x][y];


                for (size_t i = 0; i < pile->size(); i++){
                    int x = pile->at(i);
                    for (size_t j = i+1; j < pile->size(); j++){
                        int y = pile->at(j);
                        if (doublon(x, y, points, precision)){
                            a_supprimer[x] = true;
                            break;
                        }

                    }
                }
                delete pile;
            }
        }
    }

    // Suppression des doublons
    for (int i = point_counter - 1; i >= 0; i--){
        if (a_supprimer[i]){
            points->erase(points->begin() + i);
            point_counter--;
        }
    }

    // Suppression des objets
    free(a_supprimer);
    for(int x = 0; x < n; x++){
        free(grille[x]);
    }
    free(grille);
}

//Structure unir et trouver
std::vector<int>* creer_partition(int n){
    std::vector<int>* p = new std::vector<int>;
    for(int i = 0; i < n; i++){
        p->push_back(0);
    }

    for(int i = 0; i < n; i++){
        p->at(i) = i;
    }

    return p;
}

int trouver(std::vector<int>* p, int i){
    return p->at(i);
}

void unir(std::vector<int>* p, int i, int j, int n){
    int repr_i = trouver(p, i);
    int repr_j = trouver(p, j);

    for (int k = 0; k < n; k++){
        if (p->at(k) == repr_i){
            p->at(k) = repr_j;
        }
    }
}
//Fin structure unir et trouver

typedef struct {double dist;
                int i;
                int j;
            } arete;

arete* creer_arete(double d, int i, int j){
    arete* a = (arete*)malloc(sizeof(arete));
    a->dist = d;
    a->i = i;
    a->j = j;
    return a;
}

int compare_aretes(arete* a, arete* b){ return a->dist < b->dist; }

std::vector<arete*>* trier_distance(std::vector<std::shared_ptr<Point>>* points){
    int n = points->size();
    std::vector<arete*>* distances = new std::vector<arete*>;
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            double d = points->at(i)->dist(points->at(j));
            distances->push_back(creer_arete(d, i, j));
        }
    }
    std::sort(distances->begin(), distances->end(), compare_aretes);
    return distances;
}

bool est_uni(std::vector<int>* partition){
    int n = partition->size();
    for(int i = 0; i < n-1; i++){
        if (partition->at(i) != partition->at(i+1)){
            return false;
        }
    }
    return true;
}

void Quadtree::lier_points(std::vector<std::shared_ptr<Point>>* non_lies, std::vector<int>* partition, int& partition_id){
    if (divided){
        for (int i = 0; i < 4; i++){
            quadtrees[i]->lier_points(non_lies, partition, partition_id);
            partition_id++;
        }

        int n = non_lies->size();
        std::vector<arete*>* distances = trier_distance(non_lies);
        int* vus = (int*)malloc(sizeof(int)*n);
        for (int i = 0; i < n; i++){
            vus[i] = 0;
        }

        //kruskal, mais certains points sont déjà relié par la récursion
        bool flag = true;
        int indice = 0;
        while (flag && indice < n){
            arete* a = distances->at(indice);
            int i = a->i;
            int j = a->j;
            if (trouver(partition, i) !=  trouver(partition, j)){
                unir(partition, i, j, n);
                vus[i]++;
                vus[j]++;
                non_lies->at(i)->suivants->push_back(non_lies->at(j));
                non_lies->at(j)->suivants->push_back(non_lies->at(i));
            }
            if (est_uni(partition)){
                flag = false;
            }
            indice++;
        }

        for (size_t i = 0; i < distances->size(); i++){
            free(distances->at(i));
        }
        delete distances;

        std::vector<std::shared_ptr<Point>>* nouv_non_lies = new std::vector<std::shared_ptr<Point>>;
        std::vector<int>* nouv_partition = new std::vector<int>;
        partition_id++;
        for (int i = 0; i < n; i++){
            if (vus[i] == 0){
                nouv_non_lies->push_back(non_lies->at(i));
                nouv_partition->push_back(partition_id);
            }
        }
        free(vus);

        non_lies->clear();
        non_lies->insert(non_lies->end(), nouv_non_lies->begin(), nouv_non_lies->end());
        delete nouv_non_lies;

        partition->clear();
        partition->insert(partition->end(), nouv_partition->begin(), nouv_partition->end());
        delete nouv_partition;
    }
    else{
        //kruskal
        int n = points->size();

        //on veut savoir combien de fois chaque points sera lié pour l'ajouter ou non à non_lies
        int* vus = (int*)malloc(sizeof(int)*n);
        for (int i = 0; i < n; i++){
            vus[i] = 0;
        }

        std::vector<int>* p = creer_partition(n);
        std::vector<arete*>* distances = trier_distance(points);
        int compteur = 0;
        int indice = 0;
        while (compteur < n-1){
            arete* a = distances->at(indice);
            int i = a->i;
            int j = a->j;
            if (trouver(p, i) !=  trouver(p, j)){
                unir(p, i, j, n);
                vus[i]++;
                vus[j]++;
                points->at(i)->suivants->push_back(points->at(j));
                points->at(j)->suivants->push_back(points->at(i));
                compteur++;
            }
            indice++;
        }

        delete p;
        for (int i = (int)distances->size()-1 ; i >=0 ; i--){
            free(distances->at(i));
        }
        delete distances;

        partition_id++;
        for (int i = 0; i < n; i++){
            if (vus[i] == 1){
                //on indique que points[i] est dans la composante "numero_quadtree"
                non_lies->push_back(points->at(i));
                partition->push_back(partition_id);
            }
        }
        free(vus);

    }
}

int Quadtree::zone_non_vide(double x_min, double x_max, double y_min, double y_max){
    if (forme->evaluer_xy(x_min, y_min) > 0 &&
        forme->evaluer_xy(x_min, y_max) > 0 &&
        forme->evaluer_xy(x_max, y_min) > 0 &&
        forme->evaluer_xy(x_max, y_max) > 0){
            return 0;
        }
    if (forme->evaluer_xy(x_min, y_min) < 0 &&
        forme->evaluer_xy(x_min, y_max) < 0 &&
        forme->evaluer_xy(x_max, y_min) < 0 &&
        forme->evaluer_xy(x_max, y_max) < 0){
            return 0;
        }
    return 1;
}

bool Quadtree::zone_recherche(double& x_min, double& x_max, double& y_min, double& y_max){
    //On suppose que le quadtree est une feuille
    if (!divided){
        if (zone_non_vide(x_min, x_max, y_min, y_max) == 0){
            return false;
        }

        double centre_x = (x_min + x_max)/2;
        double centre_y = (y_min + y_max)/2;

        int cadran0 = zone_non_vide(x_min, centre_x, centre_y, y_max);
        int cadran1 = zone_non_vide(centre_x, x_max, centre_y, y_max);
        int cadran2 = zone_non_vide(x_min, centre_x, y_min, centre_y);
        int cadran3 = zone_non_vide(centre_x, x_max, y_min, centre_y);

        int somme = cadran0 + cadran1 + cadran2 + cadran3;

        if (somme == 3 || somme == 4){
            return true;
        }
        else if (somme == 1){
            if (cadran0 == 1){
                x_max = centre_x;
                y_min = centre_y;
                return zone_recherche(x_min, x_max, y_min, y_max);
            }
            else if (cadran1 == 1){
                x_min = centre_x;
                y_min = centre_y;
                return zone_recherche(x_min, x_max, y_min, y_max);
            }
            else if (cadran2 == 1){
                x_max = centre_x;
                y_max = centre_y;
                return zone_recherche(x_min, x_max, y_min, y_max);
            }
            else if (cadran3 == 1){
                x_min = centre_x;
                y_max = centre_y;
                return zone_recherche(x_min, x_max, y_min, y_max);
            }
        }
        else if(somme==2){
            if (cadran0 == 1){
                if (cadran1 == 1){ // haut
                    y_min = centre_y;
                    return true;
                }
                else if(cadran2 == 1){ // gauche
                    x_max = centre_x;
                    return true;
                }
            }
            else if (cadran3 == 1){
                if (cadran1 == 1){ // droite
                    x_min = centre_x;
                    return true;
                }
                else if(cadran2 == 1){
                    y_max = centre_y;
                    return true;
                }
            }
            else {
                return true;
            }
        }
    }
    return false;
    //n'arrive jamais
}

void Quadtree::copier(Quadtree* quadtree){
    //this n'est pas divisé, les quadtrees enfants n'étaient pas initialisé
    delete points;

    free(quadtrees);
    quadtrees = quadtree->quadtrees;
    points = quadtree->points;
    forme = quadtree->forme;

    min_x = quadtree->min_x;
    min_y = quadtree->min_y;
    max_x = quadtree->max_x;
    max_y = quadtree->max_y;
    center = quadtree->center;
    divided = quadtree->divided;
    
    subdivision = quadtree->subdivision;
    nb_pts_par_quadtree = quadtree->nb_pts_par_quadtree;
    point_counter = quadtree->point_counter;


    quadtree->divided = false;
    quadtree->quadtrees = nullptr;
    quadtree->points = nullptr;

    delete quadtree;
}

bool Quadtree::appartient(std::shared_ptr<Point> point){
    if (std::isnan(point->x()) || std::isnan(point->x())){
        return false;
    }
    if (point->x() > max_x || point->x() < min_x ||
        point->y() > max_y || point->y() < min_y){
        return false;
    }
    return true;
}

typedef struct triplet{ std::shared_ptr<Point> p1;
    std::shared_ptr<Point> p2;
    std::shared_ptr<Point> nouv_point;

} triplet;

std::shared_ptr<triplet> creer_triplet(std::shared_ptr<Point> p1, std::shared_ptr<Point> p2, std::shared_ptr<Point> nouv_point){
    std::shared_ptr<triplet> t = std::make_shared<triplet>();
    t->p1 = p1;
    t->p2 = p2;
    t->nouv_point = nouv_point;
    return t;
}

void lier_triplets(std::shared_ptr<Point> p1, std::shared_ptr<Point> p2, std::shared_ptr<Point> p_m){
    for(int i = p1->suivants->size()-1; i >=0 ; i--){
        if (p1->suivants->at(i).lock() == p2 || p1->suivants->at(i).expired()){
            p1->suivants->erase(p1->suivants->begin() + i);
            p1->suivants->push_back(p2);
            break;
        }
    }
    for(int i = p2->suivants->size()-1; i >=0 ; i--){
        if (p2->suivants->at(i).lock() == p1 || p2->suivants->at(i).expired()){
            p2->suivants->erase(p2->suivants->begin() + i);
            p2->suivants->push_back(p1);
            break;
        }
    }
    p_m->suivants->push_back(p1);
    p_m->suivants->push_back(p2);
}

void Quadtree::augmenter_qualite_V1_ajout(int nb_iter, double precision){
    //On suppose qu'on est arrivé dans une feuille
    if (!divided){
        std::vector<std::shared_ptr<triplet>>* pile = new std::vector<std::shared_ptr<triplet>>;
        for (size_t i = 0; i < points->size(); i++){
            std::shared_ptr<Point> point = points->at(i);

            std::shared_ptr<Point> p1 = nullptr;
            if (point->suivants->size() >= 1){
                p1 = point->suivants->at(0).lock();
            }

            std::shared_ptr<Point> p2 = nullptr;
            if (point->suivants->size() >= 2){
                p2 = point->suivants->at(1).lock();
            }  
            
            if (p1 != nullptr){
                double diff_x = std::abs(point->x() -  p1->x());
                double diff_y = std::abs(point->y() -  p1->y());

                double moy_x = (point->x() + p1->x())/2;
                double moy_y = (point->y() + p1->y())/2;

                std::shared_ptr<Point> nouv_point = std::make_shared<Point>(max_x+1, 0);
                //On l'initialise à une valeur en dehors du quadtree

                if (diff_x > diff_y){
                    Polynome* p_eval = forme->evaluer_x(moy_x);
                    Polynome* p_deriv = p_eval->deriver_y();
                    double y = p_eval->trouver_racine_newton_y(nb_iter, moy_y, p_deriv, precision);
                    nouv_point->set_x(moy_x);
                    nouv_point->set_y(y);
                    delete p_eval;
                    delete p_deriv;
                }
                else{
                    Polynome* p_eval = forme->evaluer_y(moy_y);
                    Polynome* p_deriv = p_eval->deriver_x();
                    double x = p_eval->trouver_racine_newton_x(nb_iter, moy_x, p_deriv, precision);
                    nouv_point->set_x(x);
                    nouv_point->set_y(moy_y);
                    delete p_eval;
                    delete p_deriv;
                }

                if (appartient(nouv_point)){
                    std::shared_ptr<triplet> t = creer_triplet(point, p1, nouv_point);
                    pile->push_back(t);
                }
            }

            if (p2 != nullptr){
                double diff_x = std::abs(point->x() -  p2->x());
                double diff_y = std::abs(point->y() -  p2->y());

                double moy_x = (point->x() + p2->x())/2;
                double moy_y = (point->y() + p2->y())/2;

                std::shared_ptr<Point> nouv_point = std::make_shared<Point>(max_x+1, 0);
                //On l'initialise à une valeur en dehors du quadtree

                if (diff_x > diff_y){
                    Polynome* p_eval = forme->evaluer_x(moy_x);
                    Polynome* p_deriv = p_eval->deriver_y();
                    double y = p_eval->trouver_racine_newton_y(nb_iter, moy_y, p_deriv, precision);
                    nouv_point->set_x(moy_x);
                    nouv_point->set_y(y);
                    delete p_eval;
                    delete p_deriv;
                }
                else{
                    Polynome* p_eval = forme->evaluer_y(moy_y);
                    Polynome* p_deriv = p_eval->deriver_x();
                    double x = p_eval->trouver_racine_newton_x(nb_iter, moy_x, p_deriv, precision);
                    nouv_point->set_x(x);
                    nouv_point->set_y(moy_y);
                    delete p_eval;
                    delete p_deriv;
                }

                if (appartient(nouv_point)){
                    std::shared_ptr<triplet> t = creer_triplet(point, p2, nouv_point);
                    pile->push_back(t);
                }
            }

        }

        //Ajout des nouv_points
        for (size_t i = 0; i < pile->size(); i++){
            std::shared_ptr<triplet> t = pile->at(i);
            
            lier_triplets(t->p1, t->p2, t->nouv_point);
            points->push_back(t->nouv_point);
            point_counter++;
        }
        delete pile;

        this->nettoyer_double(min_x, max_x, min_y, max_y, 10, precision);
        this->divide_space();
    }
}

void Quadtree::augmenter_qualite_V2_creation(int nb_iter, double precision, double facteur){
    if (divided){
        for (int i = 0; i < 4; i++){
            quadtrees[i]->augmenter_qualite_V2_creation(nb_iter, precision, facteur);
        }
    }
    else{
        int nb_pts = nb_pts_par_quadtree * 1 * facteur;
        // int nb_pts = point_counter * 4 * facteur;

        
        double x_min = min_x;
        double x_max = max_x;
        double y_min = min_y;
        double y_max = max_y;

        bool non_vide = zone_recherche(x_min, x_max, y_min, y_max);
        if (non_vide){
            int nb_section_doublons = 10; //arbitraire
            std::vector<std::shared_ptr<Point>>* non_lies = new std::vector<std::shared_ptr<Point>>;
            Quadtree* q = forme->generer_quadtree_poly(x_max, x_min, y_max, y_min,
                nb_iter, nb_pts, nb_pts_par_quadtree, nb_section_doublons, precision,
                forme, non_lies);
            delete non_lies;
            copier(q);
        }
        
    }
}

quad_courb* quad_courb_creer(Quadtree* q, double c){
    quad_courb* qc = (quad_courb*)malloc(sizeof(quad_courb));
    qc->quad = q;
    qc->courb = c;
    return qc;
}


void Quadtree::augmenter_qualite_visible(Camera* c, int nb_iter, double precision, int version,
                                            Quadtree* self, bool qualite, std::vector<quad_courb*>* tableau_quad_courb){
    double window_left_pos = c->camera_pos_to_world_pos_x(0);
    double window_up_pos = c->camera_pos_to_world_pos_y(0);
    double window_right_pos = c->camera_pos_to_world_pos_x(c->get_width());
    double window_down_pos = c->camera_pos_to_world_pos_y(c->get_height());

    if (divided){
        if (window_left_pos >= center->x()){
            if (window_up_pos >= center->y()){
                quadtrees[1]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[1], qualite, tableau_quad_courb);
            }
            else if (window_down_pos <= center->y()){
                quadtrees[3]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[3], qualite, tableau_quad_courb);
            }
            else {
                quadtrees[1]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[1], qualite, tableau_quad_courb);
                quadtrees[3]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[3], qualite, tableau_quad_courb);
            }
        }
        else if (window_right_pos <= center->x()){
            if (window_up_pos >= center->y()){
                quadtrees[0]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[0], qualite, tableau_quad_courb);
            }
            else if (window_down_pos <= center->y()){
                quadtrees[2]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[2], qualite, tableau_quad_courb);
            }
            else {
                quadtrees[0]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[0], qualite, tableau_quad_courb);
                quadtrees[2]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[2], qualite, tableau_quad_courb);
            }
        }
        else{
            quadtrees[0]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[0], qualite, tableau_quad_courb);
            quadtrees[1]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[1], qualite, tableau_quad_courb);
            quadtrees[2]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[2], qualite, tableau_quad_courb);
            quadtrees[3]->augmenter_qualite_visible(c, nb_iter, precision, version, quadtrees[3], qualite, tableau_quad_courb);
        }
    }
    else{
        if (version == 1){
            this->augmenter_qualite_V1_ajout(nb_iter, precision);

        }
        else if (version == 2){
            if (qualite){
                tableau_quad_courb->push_back(quad_courb_creer(self, courbure_moy_quadtree()));
            }
            else{
                this->augmenter_qualite_V2_creation(nb_iter, precision, 1);
            }
        }
    }
}

double courb(std::shared_ptr<Point> p1, std::shared_ptr<Point> p, std::shared_ptr<Point> p2){
    //p est le point au milieu
    double x_v1 = p1->x() - p->x();
    double y_v1 = p1->y() - p->y();

    double x_v2 = p2->x() - p->x();
    double y_v2 = p2->y() - p->y();

    double cross = x_v1 * y_v2 - y_v1 * x_v2; // produit vectoriel en 2D
    double dot = x_v1 * x_v2 + y_v1 * y_v2;

    return std::atan2(std::abs(cross), dot); // angle non orienté
}


double Quadtree::courbure_moy_quadtree(){
    int n = points->size();
    if (n > 0){
        double courbure = 0;
        for (int i = 0; i < n; i++){
            auto point = points->at(i);
            if (point->suivants->size() == 2){
                if(!point->suivants->at(0).expired() && !point->suivants->at(1).expired()){
                    courbure += courb(point->suivants->at(0).lock(), point, point->suivants->at(1).lock());
                }
            }  
        }
        return courbure/points->size();
    }
    return 0;
}

