#ifndef POSITION_H
#define POSITION_H
#include "matrix.h"


class Position : public Matrix
{
public:
    Position();
    Position(int Np, int Di, double b, double step);
    ~Position(){}
    void setPosition_initial();
    void setPosition(double value);
    void setRandomPosition();
    void setNewToOld(Matrix r);
    void set_step(double step){ sd = step;}
    double get_step(){return sd;}
    double get_beta(){return beta;}
    double getrsqrt();
    Matrix& get_Matrix(){return r;}
    double getPosition(int Np, int Di){return r.get_Elem(Np, Di);}
protected:
    Matrix r;
    double beta;
    double sd;




};

ostream& operator<<(ostream& os, const Matrix & r);

#endif // POSITION_H
