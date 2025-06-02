#ifndef CAM
#define CAM

#include <SFML/Graphics.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>

#include "point.h"

class Quadtree;

class Camera{
    public:
        Camera(GLFWwindow* window, int window_width, int window_height);

        void set_XY(double x, double y);
        void set_dragXY(double x, double y);
        void set_camera_dragXY(double x, double y);
        void set_afficher_points(bool val);
        void set_afficher_quadtrees(bool val);

        std::shared_ptr<Point> get_XY();
        std::shared_ptr<Point> get_dragXY();
        std::shared_ptr<Point> get_cameraXY_drag();
        bool get_flag_drag();
        double get_zoom();
        int get_mult();
        int get_width();
        int get_height();
        bool get_afficher_points();
        bool get_afficher_quadtrees();

        double camera_pos_to_world_pos_x(double x);
        double camera_pos_to_world_pos_y(double y);
        double world_pos_to_camera_pos_x(double x);
        double world_pos_to_camera_pos_y(double y);

        static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
        static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
        static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
        static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

        int render_points(Quadtree* quadtree, GLFWwindow* window);
        int render_visible_points(Quadtree* quadtree, GLFWwindow* window, Point* last_point, bool* has_last_point);
        
    private:
        double zoom;

        std::shared_ptr<Point> cameraXY;
        std::shared_ptr<Point> cameraXY_drag;
        std::shared_ptr<Point> dragXY;

        bool flag_drag;

        int MULT;
        int WINDOW_WIDTH;
        int WINDOW_HEIGHT;

        bool afficher_points;
        bool afficher_quadtrees;
};

#endif
