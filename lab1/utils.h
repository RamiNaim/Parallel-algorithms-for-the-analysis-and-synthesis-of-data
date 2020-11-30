//
// Created by Рами Ал-Наим on 30.11.2020.
//

#ifndef LAB1_BASEFUNCTION_H
#define LAB1_BASEFUNCTION_H

#endif //LAB1_BASEFUNCTION_H

class BaseFunction{
public:
    virtual double call(double x) = 0;
    virtual double integral(double a, double b) = 0;
};

template <typename T>
class Integrator{
public:
    virtual int integrate(double a, double b, T func) = 0;
};