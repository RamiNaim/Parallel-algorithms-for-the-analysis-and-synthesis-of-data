#include <iostream>
#include <omp.h>
#include "utils.h"
#include <fstream>

#define eps 0.00001

class SineFunction : public BaseFunction{
public:
    double call(double x){
        double result = 1 / (pow(x, 2.0)) * pow( sin(1/x), 2.0 );
        return result;
    }

    double integral(double a, double b){
        double result;
        result = 2 * (b - a) / (a*b) + sin(2/b) - sin(2/a);
        result *= 0.25;
        return result;
    }
};

class SerialIntegrator : public Integrator<SineFunction>{
public:
    std::string type = "Serial";
    int integrate(double a, double b, SineFunction func){
        double integral = 0, prev_integral = 0;
        double error = 1000;
        int n = 3;

        while (abs(error) > eps*integral){
            prev_integral = integral;
            integral = 0;
            double step = (b - a) / n;
            double x;
            for (int i=1; i < n-1; i++){
                x = a + i * step;
                integral += func.call(x);
            }
            integral += func.call(a)/2 + func.call(b)/2;
            integral *= step;

            error = integral - prev_integral;
            n++;
        }
        return n;
    }
};

class AtomicIntegrator : public Integrator<SineFunction>{
private:
    int numOfThreads;
public:
    std::string type = "Atomic";
    AtomicIntegrator(int numberOfThreads){
        numOfThreads = numberOfThreads;
        type += "_" + std::to_string(numOfThreads);
    }

    int getNumOfThreads(){
        return numOfThreads;
    }


    int integrate(double a, double b, SineFunction func){
        double integral = 0, prev_integral = 0;
        double error = 1000;
        int n = 3;

        while (abs(error) > eps*integral){
            prev_integral = integral;
            integral = 0;
            double step = (b - a) / n;
            double x;
            #pragma omp parallel for num_threads(numOfThreads) private(x)
            for (int i=1; i < n-1; i++){
                x = a + i * step;
                #pragma omp atomic
                integral += func.call(x);
            }
            integral += func.call(a)/2 + func.call(b)/2;
            integral *= step;

            error = integral - prev_integral;
            n++;
        }
        return n;
    }
};

template<typename T>
void run_integrals(T integrator)
{
    std::ofstream log_file("../" + integrator.type + ".csv");
    auto function = SineFunction();

    double a=0.000001, b;
    double t1,t2;
    log_file << "a,b,N,t";
    for (auto k=0;k<7;k++){
        a = a*10;
        b = a*10;

        t1 = omp_get_wtime();
        auto n_of_steps = integrator.integrate(a, b, function);
        t2 = omp_get_wtime();

        log_file << "\n" << a << "," << b << "," << n_of_steps << "," << t2 - t1;
    }

    log_file.close();

};

int main() {
    auto serial_integrator = SerialIntegrator();
    run_integrals(serial_integrator);

    auto atomic_integrator = AtomicIntegrator(6);
    run_integrals(atomic_integrator);

    return 0;
}