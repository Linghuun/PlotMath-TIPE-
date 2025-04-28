#include "camera.h"

const int MULT = 20;

Camera::Camera(GLFWwindow* window)
{
    zoom = 1.0f;
    cameraXY = new Point(0, 0);
    cameraXY_drag = new Point(0, 0);
    dragXY = new Point(0, 0);
    flag_drag = false;

    glfwSetWindowUserPointer(window, this);

    glfwSetKeyCallback(window, key_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
}

Camera::~Camera(){
    delete cameraXY;
    delete cameraXY_drag;
    delete dragXY;
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

void Camera::set_flag_drag(){ flag_drag = true; }

Point* Camera::get_XY(){ return cameraXY; }

Point* Camera::get_dragXY(){ return dragXY; }

Point* Camera::get_cameraXY_drag(){ return cameraXY_drag; }

bool Camera::get_flag_drag(){ return flag_drag; }

double Camera::get_zoom(){ return zoom; }



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

            // std::cout << "drag : " << camera->dragXY->x() << ", " << camera->dragXY->y() << std::endl;
            // std::cout << "world: " << camera->camera_pos_to_world_pos_x(camera->dragXY->x()) << ", " << camera->camera_pos_to_world_pos_y(camera->dragXY->y()) << std::endl;
            // std::cout << "cam  : " << camera->cameraXY->x() << ", " << camera->cameraXY->y() << std::endl;
            
        
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
        // std::cout << "cdrag: " << camera->cameraXY_drag->x() << ", "<<  camera->cameraXY_drag->y() << '\n';
        camera->cameraXY->set_x(camera->cameraXY_drag->x() + ((camera->dragXY->x() - xpos) / camera->zoom)/MULT);
        camera->cameraXY->set_y(camera->cameraXY_drag->y() + ((camera->dragXY->y()- ypos) / camera->zoom)/MULT);
    }
}

void Camera::scroll_callback(GLFWwindow* window, double /*xoffset*/, double yoffset) {
    Camera* camera = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    if (camera) {
        const float zoomSpeed = 0.05f;
        camera->zoom *= (1.0f + yoffset * zoomSpeed);
    }
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
    }
}


