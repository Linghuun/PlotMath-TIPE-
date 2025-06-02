#ifndef POINT
#define POINT

#include <vector>
#include <memory>

class Point : public std::enable_shared_from_this<Point> {
    public:
        Point(double x, double y);
        Point();
        ~Point();
        double dist(std::shared_ptr<Point> point);
        
        double x();
        double y();

        void set_x(double x);
        void set_y(double y);

        std::vector<std::weak_ptr<Point>>* suivants;

    protected:
        double coord_x;
        double coord_y;             
};

#endif
