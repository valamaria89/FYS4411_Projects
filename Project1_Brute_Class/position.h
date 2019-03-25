#ifndef POSITION_H
#define POSITION_H
#include "matrix.h"

class Position{

public:

    void setPosition_initial(Matrix r);
    void setPosition_(Matrix r1, Matrix r);
    void setRandomPosition(Matrix r, double step, double random);
    void setRandomStep(Matrix rold, Matrix rnew, double step, double random);
    void setNewToOld(Matrix rold, Matrix rnew);
    void set_h(Matrix r, double h);
    double getrsqrt(Matrix r);



};



#endif // POSITION_H
