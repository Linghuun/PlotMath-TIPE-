#include "camera.h"
#include "quadtree.h"
#include "polynome.h"

#include <SFML/Graphics.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <memory>


const int NB_MINI_PTS = 100;
const int NB_MAX_PTS = 500;

const double FACTEUR_MIN = 0.75;
const double FACTEUR_MAX = 2;


double courb_moy(std::vector<quad_courb*>* tableau_quad_courb){
    double courbure = 0;
    int nb_pts = 0;
    for (size_t i = 0; i < tableau_quad_courb->size(); i++){
        quad_courb* qc = tableau_quad_courb->at(i);
        courbure += qc->courb * qc->quad->point_counter;
        nb_pts += qc->quad->point_counter;
    }
    if (nb_pts > 0){
        return courbure/nb_pts;
    }

    printf("Problème courb_moy()\n");
    return NAN;
}

double courb_max(std::vector<quad_courb*>* tableau_quad_courb){
    double courbure_max = 0;
    for (size_t i = 0; i < tableau_quad_courb->size(); i++){
        quad_courb* qc = tableau_quad_courb->at(i);
        if (qc->courb > courbure_max){
            courbure_max = qc->courb;
        }
    }
    if (courbure_max > 0){
        return courbure_max;
    }

    printf("Problème courb_max()\n");
    return NAN;
}

double facteur(double courbure_moyenne, double courbure, double courbure_max){
    if (courbure_moyenne == courbure_max){
        return 1;
    }
    if (courbure < courbure_moyenne){
        return FACTEUR_MIN + courbure * ((1-FACTEUR_MIN)/(courbure_moyenne));
    }
    //else
    return 1 + (courbure - courbure_moyenne) * ((FACTEUR_MAX-1)/(courbure_max-courbure_moyenne));
}

int main(int /*argc*/, char *argv[]) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, 255); // define the range

    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    int WINDOW_WIDTH = 1000;
    int WINDOW_HEIGHT = 960;
    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "My Title", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window!" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW!" << std::endl;
        return -1;
    }

    // Set the viewport to match the window dimensions
    glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    // Set up the projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, WINDOW_WIDTH, 0.0, WINDOW_HEIGHT, -1.0, 1.0);

    // Set up the model-view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Create the camera assigned to the current window
    Camera* camera = new Camera(window, WINDOW_WIDTH, WINDOW_HEIGHT);

    Polynome* forme = new Polynome();

    //Polynome lisse
    forme->ajouter_terme(std::unique_ptr<Terme>(new Terme(6, 0, 2)));
    forme->ajouter_terme(std::unique_ptr<Terme>(new Terme(0, 6, 1)));
    forme->ajouter_terme(std::unique_ptr<Terme>(new Terme(0, 2, 2)));
    forme->ajouter_terme(std::unique_ptr<Terme>(new Terme(0, 3, 1)));
    forme->ajouter_terme(std::unique_ptr<Terme>(new Terme(3, 0, 4)));
    forme->ajouter_terme(std::unique_ptr<Terme>(new Terme(2, 1, -3)));
    forme->ajouter_terme(std::unique_ptr<Terme>(new Terme(0, 2, -3)));    

    double precision = 10e-7;
    int nb_pts = std::stoi(argv[3]);
    int nb_iter = 5;
    int nb_pts_par_quadtree = std::stoi(argv[4]);
    int nb_section_doublons = 1 + (int)nb_pts*(3/4);

    auto non_lies = new std::vector<std::shared_ptr<Point>>;
    Quadtree* quadtree = forme->generer_quadtree_poly(1.5, -1.5, 1.5, -1.5, nb_iter, nb_pts, nb_pts_par_quadtree, nb_section_doublons, precision, forme, non_lies);

    if (non_lies->size() == 2){
        non_lies->at(0)->suivants->push_back(non_lies->at(1));
        non_lies->at(0)->suivants->push_back(non_lies->at(1));
    }
    std::cout << "nb_points :" << quadtree->point_counter << "\n";
    delete non_lies;

    // Boucle principale

    //render_visible_points a besoin d'un point de départ inutile
    bool has_last_point = false;
    Point* last_point = new Point(0,0);
    int nb_pts_visible = 0;
    int nb_frames = 0;
    while (!glfwWindowShouldClose(window)) {
        nb_frames++;
        if (nb_frames % 60 == 0){ //on affiche le nombre de points toutes les 2 secondes
            std::cout << "nb_pts_visible : " << nb_pts_visible << "\n";
        }

        // Nettoyer l'écran
        glClear(GL_COLOR_BUFFER_BIT);

        // Rendre les points
        nb_pts_visible = camera->render_visible_points(quadtree, window, last_point, &has_last_point);

        // augmenter_qualite
        if (nb_pts_visible < NB_MINI_PTS){
            if (std::stoi(argv[1]) == 1){
                quadtree->augmenter_qualite_visible(camera, nb_iter+10, precision*0.1, std::stoi(argv[1]), nullptr, std::stoi(argv[2]), nullptr); 
            }
            else if (std::stoi(argv[1]) == 2){
                if (std::stoi(argv[2]) == 1){
                    auto tableau_quad_courb = new std::vector<quad_courb*>;
                    quadtree->augmenter_qualite_visible(camera, nb_iter+10, precision*0.1, std::stoi(argv[1]), quadtree, std::stoi(argv[2]), tableau_quad_courb); 
                    double courbure_moyenne = courb_moy(tableau_quad_courb);
                    double courbure_max = courb_max(tableau_quad_courb);
                    for (size_t i = 0; i < tableau_quad_courb->size(); i++){
                        quad_courb* qc = tableau_quad_courb->at(i);
                        qc->quad->augmenter_qualite_V2_creation(nb_iter+10, precision*0.1, facteur(courbure_moyenne, qc->courb, courbure_max));
                    }

                    for (size_t i = 0; i < tableau_quad_courb->size(); i++){
                        free(tableau_quad_courb->at(i));
                    }
                    delete tableau_quad_courb;
                }
                else{
                    quadtree->augmenter_qualite_visible(camera, nb_iter+10, precision*0.1, std::stoi(argv[1]), quadtree, std::stoi(argv[2]), nullptr); 
                }
            }
        }
        
        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }
    std::cout << "nb_frames : " << nb_frames << "\n";
    delete last_point;

    // Clean up and close
    glfwDestroyWindow(window);
    glfwTerminate();

    delete camera;
    delete quadtree;
    delete forme;
    return 0;
}
