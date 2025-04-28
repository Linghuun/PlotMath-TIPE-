#ifndef POINT
#define POINT

class Point{
    public:
        Point(double x, double y);
        Point();
        double dist(Point* point);
        
        double x();
        double y();

        void set_x(double x);
        void set_y(double y);

        Point* previous;
        Point* next;

    protected:
        double coord_x;
        double coord_y;             
};

#endif