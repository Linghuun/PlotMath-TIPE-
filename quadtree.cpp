#include "quadtree.h"
#include "polynome.h"

#include <cstdlib>
#include <iostream>
#include <float.h>
#include <algorithm>
#include <cmath>

Quadtree::Quadtree(int subdiv, int nb_pts_max){
    quadtrees = (Quadtree**)malloc(sizeof(Quadtree*)*subdiv);

    divided = false;
    max_x = -DBL_MAX;
    max_y = -DBL_MAX;
    min_x = DBL_MAX;
    min_y = DBL_MAX;
    center = new Point(0,0);
    point_counter = 0;
    subdivision = subdiv;
    nb_pts_par_quadtree = nb_pts_max;
    points = new std::vector<Point*>;
    courbure_init = false;
}

Quadtree::~Quadtree(){
    if (divided){
        for (int i = 0; i < subdivision; i++){
            delete quadtrees[i];
        }
    }
    else{
        for (size_t i = 0; i < points->size(); i++){
            delete (*points)[i];
        }
    }
    free(quadtrees);
    delete points;
    delete center;

    if (courbure_init){
        for (int i = 0; i < 5; i++){
            delete courbure[i];
        }
        free(courbure);
    }
}

void Quadtree::add_point(Point* point){
    if(point->x()<min_x){min_x = point->x();}
    if(point->x()>max_x){max_x = point->x();}
    if(point->y()<min_y){min_y = point->y();}
    if(point->y()>max_y){max_y = point->y();}
    
    
    point_counter++;
    points->push_back(point);

}

void Quadtree::push_vector(std::vector<Point*>* points){
    for (size_t i = 0; i < points->size(); i++){
        add_point((*points)[i]);
    }
}

int Quadtree::get_point_count(){
    return point_counter;
}

Point* Quadtree::get_point(int i){
    return (*points)[i];
}

bool Quadtree::is_divided(){
    return divided;
}

Point* Quadtree::get_center(){
    return center;
}

void Quadtree::divide_space(){
    if (point_counter >= nb_pts_par_quadtree){
        center->set_x((max_x+min_x)/2);
        center->set_y((max_y+min_y)/2),
        divided = true;
        for (int i = 0; i < subdivision; i++){
            quadtrees[i] = new Quadtree(subdivision, nb_pts_par_quadtree);
        }
        quadtrees[0]->min_x = min_x;
        quadtrees[0]->max_y = max_y;
        quadtrees[0]->max_x = center->x();
        quadtrees[0]->min_y = center->y();

        quadtrees[1]->min_x = center->x();
        quadtrees[1]->max_y = max_y;
        quadtrees[1]->max_x = max_x;
        quadtrees[1]->min_y = center->y();

        quadtrees[2]->min_x = min_x;
        quadtrees[2]->max_y = center->y();
        quadtrees[2]->max_x = center->x();
        quadtrees[2]->min_y = min_y;

        quadtrees[3]->min_x = center->x();
        quadtrees[3]->max_y = center->y();
        quadtrees[3]->max_x = max_x;
        quadtrees[3]->min_y = min_y;

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
        quadtrees[0]->divide_space();
        quadtrees[1]->divide_space();
        quadtrees[2]->divide_space();
        quadtrees[3]->divide_space();
    }
}

int_pile* init_pile(int val){
    int_pile* pile = (int_pile*)malloc(sizeof(int_pile));
    pile->val = val;
    pile->next = NULL;
    return pile;
}

void liberer_int_pile(int_pile* pile){
    if (pile->next != NULL){
        liberer_int_pile(pile->next);
    }
    free(pile);
}

int_pile* ajouter_val(int_pile* pile, int val){
    int_pile* new_pile = init_pile(val);
    new_pile->next = pile;
    return new_pile;
}

void Quadtree::nettoyer_double(double precision){
    //On veut faire une grille de hachage de section espilon
    
    int_pile* pile = NULL;
    for (int i = 0; i < point_counter; i++){
        for (int j = i+1; j < point_counter; j++){
            if ((*points)[i]->dist((*points)[j]) < precision){
                pile = ajouter_val(pile, i);
                break;
            }
        }
    }
    int_pile* pointeur = pile;
    while (pointeur != NULL){
        int indice = pointeur->val;
        delete (*points)[indice];
        points->erase(points->begin() + indice); //transfome l'indice en iterator
        pointeur = pointeur->next;
    }
    printf("\n");
    liberer_int_pile(pile);
}

double plus_proche_voisin(Point* point, std::vector<Point*> *points, bool* vus, int n, int* indice){
    double distance_min = DBL_MAX;
    for (int i = 0; i < n; i++){
        if (!vus[i]){
            double distance = point->dist((*points)[i]);
            if (distance < distance_min){
                distance_min = distance;
                *indice = i;
            }
        }
    }
    return distance_min;
}

void lier_points_rec(std::vector<Point*> *points, bool* vus, int indice1, int indice2, int n, std::vector<Point*>* non_lies){
    vus[indice1] = true;
    vus[indice2] = true;

    //création des deux candidats (un seul sera retenu)
    int candidat1 = -1;
    double dist1 = plus_proche_voisin((*points)[indice1], points, vus, n, &candidat1);

    int candidat2 = -1;
    double dist2 = plus_proche_voisin((*points)[indice2], points, vus, n, &candidat2);

    if (candidat1 != -1){ // si on passe, on a aussi candidat2 != -1
        if (dist1 < dist2){

            (*points)[indice1]->previous = (*points)[candidat1];
            (*points)[candidat1]->next = (*points)[indice1];
            lier_points_rec(points, vus, candidat1, indice2, n, non_lies);
        }
        else{
            (*points)[candidat2]->previous = (*points)[ indice2];
            (*points)[ indice2]->next = (*points)[candidat2];
            lier_points_rec(points, vus, indice1, candidat2, n, non_lies);
        }
    }
    else{ //alors on a aussi candidat2 = -1 tous les points ont été vus   
        non_lies->push_back((*points)[indice1]);
        non_lies->push_back((*points)[indice2]);
    }
}

bool compare (int i,int j) { return (i>j); }

void lier_non_lies(std::vector<Point*>* non_lies, double pas_max){
    int taille = non_lies->size();
    if (taille > 1){
        std::vector<int>* val_a_suppr = new std::vector<int>;
        for (int i = 0; i < taille; i++){
            double dist_min = DBL_MAX;
            int indice = -1;
            for (int  j = 0; j < taille; j++){
                if (i != j){
                    double dist = (*non_lies)[i]->dist((*non_lies)[j]);
                    if (dist < dist_min){
                        dist_min = dist;
                        indice = j;
                    }
                }
            }
            if (dist_min < pas_max){
                
                if ((*non_lies)[i]->previous == NULL){
                    (*non_lies)[i]->previous = (*non_lies)[indice];
                    if ((*non_lies)[indice]->previous == NULL){
                        (*non_lies)[indice]->previous = (*non_lies)[i];
                    }
                    else{
                        (*non_lies)[indice]->next = (*non_lies)[i];
                    }
                }
                else{
                    (*non_lies)[i]->next = (*non_lies)[indice];
                    if ((*non_lies)[indice]->previous == NULL){
                        (*non_lies)[indice]->previous = (*non_lies)[i];
                    }
                    else{
                        (*non_lies)[indice]->next = (*non_lies)[i];
                    }
                }
                val_a_suppr->push_back(i);
                /*On n'ajoute pas 'indice', il passera plus tard*/
            }
        }
        std::sort (val_a_suppr->begin(), val_a_suppr->end(), compare);
        for(size_t i = 0; i < val_a_suppr->size(); i++){
            non_lies->erase(non_lies->begin() + (*val_a_suppr)[i]);
        }
        delete val_a_suppr;
    }
}

void Quadtree::lier_points(std::vector<Point*>* non_lies, double pas_max){
    if (is_divided()){
        for (int i = 0; i < subdivision; i++){
            quadtrees[i]->lier_points(non_lies, pas_max);
        }
        lier_non_lies(non_lies, pas_max);
    }
    else{
        
        if (points->size()>1){
            int depart = 0; 
            // choix arbitraire, il faut juste commence quelque part

            bool* vus = (bool*)malloc(sizeof(bool)*points->size());
            for (size_t i = 0; i < points->size(); i++){
                vus[i] = false;
            }
            vus[depart] = true;

            int indice1 = -1;
            plus_proche_voisin((*points)[depart], points, vus, points->size(), &indice1);
            vus[indice1] = true;

            int indice2 = -1;
            plus_proche_voisin((*points)[depart], points, vus, points->size(), &indice2);

            if (indice2 == -1){ //alors il n'y avait que 2 points
                (*points)[indice1]->previous = (*points)[depart];
                (*points)[depart]->next = (*points)[indice1];
                non_lies->push_back((*points)[indice1]);
                non_lies->push_back((*points)[depart]);
            }
            else{
                //on lie les points au départ
                (*points)[depart]->previous = (*points)[indice1];
                (*points)[depart]->next = (*points)[indice2];
                (*points)[indice2]->previous = (*points)[depart];
                (*points)[indice1]->next = (*points)[depart];

                lier_points_rec(points, vus, indice1, indice2, points->size(), non_lies);
            }
            free(vus);
        }
        else{
            if (points->size() == 1){
                non_lies->push_back((*points)[0]);
            }
        }
    }
}

double Quadtree::calculer_courbure_xy(double x, double y){
    double py = courbure[0]->evaluer_xy(x, y);
    double px = courbure[1]->evaluer_xy(x, y);
    double pyy = courbure[2]->evaluer_xy(x, y);
    double pxx = courbure[3]->evaluer_xy(x, y);
    double pxy = courbure[4]->evaluer_xy(x, y);
    return (py*py*pxx + px*px*pyy -2*py*px*pxy)/(py*py + px*px)*std::sqrt(py*py + px*px);
}