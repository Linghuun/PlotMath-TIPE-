#include "camera.h"
#include "quadtree.h"

Camera::Camera(GLFWwindow* window, int window_width, int window_height)
{
    zoom = 1.0f;
    cameraXY = std::make_shared<Point>(0, 0);
    cameraXY_drag = std::make_shared<Point>(0, 0);
    dragXY = std::make_shared<Point>(0, 0);
    flag_drag = false;
    afficher_points = false;
    afficher_quadtrees = false;

    MULT = 200;
    WINDOW_WIDTH = window_width;
    WINDOW_HEIGHT = window_height;

    glfwSetWindowUserPointer(window, this);

    glfwSetKeyCallback(window, key_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
}

void Camera::set_XY(double x, double y){
    cameraXY->set_x(x);
    cameraXY->set_y(y);
}

void Camera::set_dragXY(double x, double y){
    dragXY->set_x(x);
    dragXY->set_y(y);
}

void Camera::set_camera_dragXY(double x, double y){
    cameraXY_drag->set_x(x);
    cameraXY_drag->set_y(y);
}

void Camera::set_afficher_points(bool val){ afficher_points = val; }

void Camera::set_afficher_quadtrees(bool val){ afficher_quadtrees = val; }

std::shared_ptr<Point> Camera::get_XY(){ return cameraXY; }
std::shared_ptr<Point> Camera::get_dragXY(){ return dragXY; }
std::shared_ptr<Point> Camera::get_cameraXY_drag(){ return cameraXY_drag; }
bool Camera::get_flag_drag(){ return flag_drag; }

double Camera::get_zoom(){ return zoom; }

int Camera::get_mult(){ return MULT; }

int Camera::get_width(){ return WINDOW_WIDTH; }

int Camera::get_height(){ return WINDOW_HEIGHT; }

bool Camera::get_afficher_points(){ return afficher_points; }

bool Camera::get_afficher_quadtrees(){ return afficher_quadtrees; }


double Camera::camera_pos_to_world_pos_x(double c_point_x){
    return (c_point_x/zoom)/MULT + cameraXY->x();
}

double Camera::camera_pos_to_world_pos_y(double c_point_y){
    return (c_point_y/zoom)/MULT - cameraXY->y(); 
}

double Camera::world_pos_to_camera_pos_x(double point_x){
    return MULT*((point_x - cameraXY->x()) * zoom);
}

double Camera::world_pos_to_camera_pos_y(double point_y){
    return MULT*((point_y + cameraXY->y()) * zoom);
}

void Camera::mouse_button_callback(GLFWwindow* window, int button, int action, int /*mods*/) {
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    if (camera) {
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {

            double dragX, dragY;
            glfwGetCursorPos(window, &dragX, &dragY);
            camera->dragXY->set_x(dragX);
            camera->dragXY->set_y(dragY);

            camera->cameraXY_drag->set_x(camera->cameraXY->x());
            camera->cameraXY_drag->set_y(camera->cameraXY->y());
            camera->flag_drag = true;
        }
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
            camera->flag_drag = false;
        }
    }
}

void Camera::cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    if (camera && camera->flag_drag == true) {
        camera->cameraXY->set_x(camera->cameraXY_drag->x() + ((camera->dragXY->x() - xpos) / camera->zoom)/camera->get_mult());
        camera->cameraXY->set_y(camera->cameraXY_drag->y() + ((camera->dragXY->y()- ypos) / camera->zoom)/camera->get_mult());
    }
}

void Camera::scroll_callback(GLFWwindow* window, double /*xoffset*/, double yoffset) {
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    if (camera) {
        const float zoomSpeed = 0.05f;
        camera->zoom *= (1.0f + yoffset * zoomSpeed);
    }
}

void actionQ() {
    std::cout << "Touche Q pressée !" << std::endl;
}

void actionP() {
    std::cout << "Touche P pressée !" << std::endl;
}

void Camera::key_callback(GLFWwindow* window, int key, int /*scancode*/, int action, int /*mods*/) {
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    if (camera) {
        const float moveSpeed = 10.0f / camera->zoom;

        if (action == GLFW_PRESS || action == GLFW_REPEAT) {
            switch (key) {
                case GLFW_KEY_LEFT:
                    camera->cameraXY->set_x(camera->cameraXY->x() - moveSpeed);
                    break;
                case GLFW_KEY_RIGHT:
                    camera->cameraXY->set_x(camera->cameraXY->x() + moveSpeed);
                    break;
                case GLFW_KEY_UP:
                    camera->cameraXY->set_y(camera->cameraXY->y() + moveSpeed);
                    break;
                case GLFW_KEY_DOWN:
                    camera->cameraXY->set_y(camera->cameraXY->y() - moveSpeed);
                    break;
            }
        }

        if (action == GLFW_PRESS){
            if (key == GLFW_KEY_A){//Touche Q car clavier inversé
                camera->set_afficher_quadtrees(! camera->get_afficher_quadtrees());
            }
            else if (key == GLFW_KEY_P){
                camera->set_afficher_points(! camera->get_afficher_points());
            }
        }
    }
}

int Camera::render_points(Quadtree* quadtree, GLFWwindow* window){
    if (quadtree->divided){
        int nb_pts = 0;
        for (int i = 0; i < 4; i++){
            nb_pts += render_points(quadtree->quadtrees[i], window);
        }
        return nb_pts;
    }
    else{
        glPointSize(8.0f);
        glLineWidth(2.0f);
        glColor3f(0.0, 0.0, 0.0);
        for (int i = 0; i < quadtree->point_counter; i++) {
            auto p = quadtree->points->at(i);
            float screenX = world_pos_to_camera_pos_x(p->x());
            float screenY = world_pos_to_camera_pos_y(p->y());

            if (afficher_points){
                glBegin(GL_POINTS);
                    glVertex2f(screenX, screenY);
                glEnd();
            }

            
            for (size_t j = 0; j < p->suivants->size(); j++){
                std::shared_ptr<Point> suivant = p->suivants->at(j).lock();
                if (suivant != nullptr){
                    float screenX_suivant = world_pos_to_camera_pos_x(suivant->x());
                    float screenY_suivant = world_pos_to_camera_pos_y(suivant->y());
                    glBegin(GL_LINES);
                        glVertex2f(screenX, screenY); 
                        glVertex2f(screenX_suivant, screenY_suivant);
                    glEnd();
                }
            }
        }

        if (afficher_quadtrees){
            glColor3f(1.0, 0.0, 0.0);
            glBegin(GL_LINES);
                glVertex2f(world_pos_to_camera_pos_x(quadtree->min_x), world_pos_to_camera_pos_y(quadtree->min_y)); 
                glVertex2f(world_pos_to_camera_pos_x(quadtree->min_x), world_pos_to_camera_pos_y(quadtree->max_y));
            glEnd();

            glBegin(GL_LINES);
                glVertex2f(world_pos_to_camera_pos_x(quadtree->max_x), world_pos_to_camera_pos_y(quadtree->min_y)); 
                glVertex2f(world_pos_to_camera_pos_x(quadtree->max_x), world_pos_to_camera_pos_y(quadtree->max_y));
            glEnd();

            glBegin(GL_LINES);
                glVertex2f(world_pos_to_camera_pos_x(quadtree->min_x), world_pos_to_camera_pos_y(quadtree->min_y)); 
                glVertex2f(world_pos_to_camera_pos_x(quadtree->max_x), world_pos_to_camera_pos_y(quadtree->min_y));
            glEnd();

            glBegin(GL_LINES);
                glVertex2f(world_pos_to_camera_pos_x(quadtree->min_x), world_pos_to_camera_pos_y(quadtree->max_y)); 
                glVertex2f(world_pos_to_camera_pos_x(quadtree->max_x), world_pos_to_camera_pos_y(quadtree->max_y));
            glEnd();
        }
        
        return quadtree->point_counter;
    }
}

int Camera::render_visible_points(Quadtree* quadtree, GLFWwindow* window, Point* last_point, bool* has_last_point){
    double window_left_pos = camera_pos_to_world_pos_x(0);
    double window_up_pos = camera_pos_to_world_pos_y(0);
    double window_right_pos = camera_pos_to_world_pos_x(get_width());
    double window_down_pos = camera_pos_to_world_pos_y(get_height());
    
    if (quadtree->divided){
        int nb_pts = 0;
        if (window_left_pos >= quadtree->center->x()){
            if (window_up_pos >= quadtree->center->y()){
                nb_pts += render_visible_points(quadtree->quadtrees[1], window, last_point, has_last_point);
            }
            else if (window_down_pos <= quadtree->center->y()){
                nb_pts += render_visible_points(quadtree->quadtrees[3], window, last_point, has_last_point);
            }
            else {
                nb_pts += render_visible_points(quadtree->quadtrees[1], window, last_point, has_last_point);
                nb_pts += render_visible_points(quadtree->quadtrees[3], window, last_point, has_last_point);
            }
        }
        else if (window_right_pos <= quadtree->center->x()){
            if (window_up_pos >= quadtree->center->y()){
                nb_pts += render_visible_points(quadtree->quadtrees[0], window, last_point, has_last_point);
            }
            else if (window_down_pos <= quadtree->center->y()){
                nb_pts += render_visible_points(quadtree->quadtrees[2], window, last_point, has_last_point);
            }
            else {
                nb_pts += render_visible_points(quadtree->quadtrees[0], window, last_point, has_last_point);
                nb_pts += render_visible_points(quadtree->quadtrees[2], window, last_point, has_last_point);
            }
        }
        else{
            nb_pts += render_visible_points(quadtree->quadtrees[0], window, last_point, has_last_point);
            nb_pts += render_visible_points(quadtree->quadtrees[1], window, last_point, has_last_point);
            nb_pts += render_visible_points(quadtree->quadtrees[2], window, last_point, has_last_point);
            nb_pts += render_visible_points(quadtree->quadtrees[3], window, last_point, has_last_point);
        }
        return nb_pts;
    }
    else{
        return render_points(quadtree, window);
    }
}
