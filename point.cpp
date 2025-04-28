#include "point.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

Point::Point(double x, double y){
    coord_x = x;
    coord_y = y;
    previous = NULL;
    next = NULL;
}

Point::Point(){
    coord_x = 0;
    coord_y = 0;
    previous = NULL;
    next = NULL;
}

double Point::dist(Point* point){
    if (point == NULL) {return -1;}
    return sqrt(std::pow(point->x()-coord_x, 2) + std::pow(point->y()-coord_y, 2));
}

double Point::x(){ return coord_x; }

double Point::y(){ return coord_y; }

void Point::set_x(double x_){ coord_x = x_; }

void Point::set_y(double y_){ coord_y = y_; }