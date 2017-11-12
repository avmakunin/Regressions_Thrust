#ifndef REGRESSIONS_THRUST_KOEFFS_H
#define REGRESSIONS_THRUST_KOEFFS_H

#include <thrust/device_vector.h>

#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <fstream>
#include <ctime>
#include <sstream>
using namespace std;


struct num
{
    double c1,c2;
    num(double c1, double c2) : c1(c1), c2(c2){};

    __host__ __device__
    double operator()(double x, double y)
    {
        return (x-c1)*(y-c2);
    }
};

struct den
{
    double c1;
    den(double c1) : c1(c1){};

    __host__ __device__
    double operator()(double x)
    {
        return (x-c1)*(x-c1);
    }
};


struct findlog
{
    double Ai,Bi,Di;
    findlog(double Bi, double Di) : Bi(Bi), Di(Di) {};
    __host__ __device__
    double operator()(double x)
    {
        return log(1 + Bi*Di*x);
    }

};
struct findS
{
    double Ai,Bi;
    findS(double Ai, double Bi) : Ai(Ai), Bi(Bi){};

    __host__ __device__
    double operator()(double x, double y)
    {
        return (x - Ai - (1. / Bi) * y)*(x - Ai - (1. / Bi) * y);
    }
};


struct findS4
{
    double Y,Ai,Bi,Di,tau,a;
    findS4(double Ai, double Bi, double Di, double tau, double a) : Ai(Ai), Bi(Bi), Di(Di), tau(tau), a(a) {};

    __host__ __device__
    double operator()(double x, double y)
    {
        if (x <= tau) Y = Ai - (1. / Bi) * log(1 + x);
        else Y = Ai - (1. / Bi)*log(1 + Bi*a*(x - tau)*pow(1 + tau, -a - 1));
        return (y - Y)*(y - Y);
    }
};

#endif //REGRESSIONS_THRUST_KOEFFS_H
