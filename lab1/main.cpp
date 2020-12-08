#include <iostream>
#include <omp.h>
#include "utils.h"
#include <fstream>

#define N 100000

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
private:
    int numOfThreads;
public:
    std::string type = "Serial";
    SerialIntegrator(){
        numOfThreads = 1;
    }

    int getNumOfThreads(){
        return numOfThreads;
    }

    int integrate(double a, double b, SineFunction func, int n_steps = N){
        double integral = 0;

        double step = (b - a) / n_steps;
        double x;
        for (int i=1; i < n_steps; i++){
            x = a + i * step;
            integral += func.call(x);
        }
        integral += func.call(a)/2 + func.call(b)/2;
        integral *= step;

        double true_integral = func.integral(a, b);
        return abs(integral - true_integral) / true_integral ;
    }
};

class AtomicIntegrator : public Integrator<SineFunction>{
private:
    int numOfThreads;
public:
    std::string type = "Atomic";
    AtomicIntegrator(int numberOfThreads){
        numOfThreads = numberOfThreads;
        // type += "_" + std::to_string(numOfThreads);
    }

    int getNumOfThreads(){
        return numOfThreads;
    }


    int integrate(double a, double b, SineFunction func, int n_steps = N){
        double integral = 0;

        double step = (b - a) / n_steps;
        double x;
        #pragma omp parallel for num_threads(numOfThreads) private(x)
        for (int i=0; i < n_steps; i++){
            x = a + i * step;
            #pragma omp atomic
            integral += func.call(x);
        }
        integral += func.call(a)/2 + func.call(b)/2;
        integral *= step;

        double true_integral = func.integral(a, b);
        return abs(integral - true_integral) / true_integral ;
    }
};

class CriticalIntegrator : public Integrator<SineFunction>{
private:
    int numOfThreads;
public:
    std::string type = "Critical";

    CriticalIntegrator(int numberOfThreads){
        numOfThreads = numberOfThreads;
        // type += "_" + std::to_string(numOfThreads);
    }

    int getNumOfThreads(){
        return numOfThreads;
    }


    int integrate(double a, double b, SineFunction func, int n_steps=N){
        double integral = 0;

        double step = (b - a) / n_steps;
        double x;
        #pragma omp parallel for num_threads(numOfThreads) private(x)
        for (int i=0; i < n_steps; i++){
            x = a + i * step;
            #pragma omp critical{
                integral += func.call(x);
            }
        }
        integral += func.call(a)/2 + func.call(b)/2;
        integral *= step;

        double true_integral = func.integral(a, b);
        return abs(integral - true_integral) / true_integral ;
    }
};

class LockIntegrator : public Integrator<SineFunction>{
private:
    int numOfThreads;
public:
    std::string type = "Lock";

    LockIntegrator(int numberOfThreads){
        numOfThreads = numberOfThreads;
        // type += "_" + std::to_string(numOfThreads);
    }

    int getNumOfThreads(){
        return numOfThreads;
    }

    int integrate(double a, double b, SineFunction func, int n_steps=N){
        double integral = 0;

        omp_lock_t integral_update_lock;
        omp_init_lock(&integral_update_lock);
        integral = 0;
        double step = (b - a) / n_steps;
        double x;
        #pragma omp parallel for num_threads(numOfThreads) private(x)
        for (int i=0; i < n_steps; i++){
            x = a + i * step;
            omp_set_lock(&integral_update_lock);
            integral += func.call(x);
            omp_unset_lock(&integral_update_lock);
        }
        integral += func.call(a)/2 + func.call(b)/2;
        integral *= step;

        double true_integral = func.integral(a, b);
        return abs(integral - true_integral) / true_integral ;
    }
};

class ReductionIntegrator : public Integrator<SineFunction>{
private:
    int numOfThreads;
public:
    std::string type = "Reduction";

    ReductionIntegrator(int numberOfThreads){
        numOfThreads = numberOfThreads;
        // type += "_" + std::to_string(numOfThreads);
        // omp_set_num_threads(numOfThreads);
    }

    int getNumOfThreads(){
        return numOfThreads;
    }

    int integrate(double a, double b, SineFunction func, int n_steps = N){
        double integral = 0;

        double step = (b - a) / n_steps;
        double x;
        #pragma omp parallel for private(x) reduction(+:integral)
        for (int i = 1; i < n_steps; i++) {
            x = a + i * step;
            integral += func.call(x);
        }
        integral += func.call(a)/2 + func.call(b)/2;
        integral *= step;

        double true_integral = func.integral(a, b);
        return abs(integral - true_integral) / true_integral ;
    }
};

template<typename T>
void run_integrals(T integrator){

    std::ofstream log_file("../result/" + integrator.type + "/" + integrator.type + "_" +
                           std::to_string(integrator.getNumOfThreads()) + ".csv");

    auto function = SineFunction();

    std::cout << "Start " << integrator.type << " computation";
    if (integrator.getNumOfThreads() != 1){
        std::cout << " with " << integrator.getNumOfThreads() << " treads" << std::endl;
    }
    else{
        std::cout << std::endl;
    }

    double a=0.000001, b;
    double t1,t2;
    log_file << "a,b,eps,t";
    for (auto k=0;k<7;k++){
        a = a*10;
        b = a*10;

        t1 = omp_get_wtime();
        auto eps = integrator.integrate(a, b, function);
        t2 = omp_get_wtime();

        log_file << "\n" << a << "," << b << "," << eps << "," << t2 - t1;
    }

    log_file.close();

};

int main() {

    auto serial_integrator = SerialIntegrator();
    run_integrals(serial_integrator);

    for (auto i=1; i<=8; i++) {

        int n_threads = i*2;
        auto atomic_integrator = AtomicIntegrator(n_threads);
        run_integrals(atomic_integrator);

        auto critical_integrator = CriticalIntegrator(n_threads);
        run_integrals(critical_integrator);

        auto lock_integrator = LockIntegrator(n_threads);
        run_integrals(lock_integrator);

        auto reduction_integrator = ReductionIntegrator(n_threads);
        run_integrals(reduction_integrator);
    }

    return 0;
}