#ifndef CAM
#define CAM

#include "quadtree.h"

#include <SFML/Graphics.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>

class Camera{
    public:
        Camera(GLFWwindow* window);
        ~Camera();

        void set_XY(double x, double y);
        void set_dragXY(double x, double y);
        void set_camera_dragXY(double x, double y);
        void set_flag_drag();

        Point* get_XY();
        Point* get_dragXY();
        Point* get_cameraXY_drag();
        bool get_flag_drag();
        double get_zoom();

        double camera_pos_to_world_pos_x(double x);
        double camera_pos_to_world_pos_y(double y);
        double world_pos_to_camera_pos_x(double x);
        double world_pos_to_camera_pos_y(double y);

    
        static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
        static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
        static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
        static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);


    private:
        double zoom;

        Point* cameraXY;

        Point* cameraXY_drag;

        Point* dragXY;

        bool flag_drag;
};

#endif