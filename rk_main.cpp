#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include "runge_kutta.hpp"

using double_array = std::array<double, 2>;

double dy_dx(double x, const double_array& a)
{
    (void)x;
    return a[1];
}

double dz_dx(double x, const double_array& a)
{
    return (cos(3.0 * x) - 4.0 * a[0]);
}

int main(int argc, char** argv)
{
    std::cout << "Runge-Kutta method\n";

    double_array inits = {0.8, 2.};
    std::array<std::function<double(double, const double_array&)>, 2> differential_equations = {dy_dx, dz_dx};

    // vectors for result storing
    //    std::vector<double> y{};
    //    std::vector<double> y_first_derivative{};
    //    std::array<std::back_insert_iterator<std::vector<double>>, 2> output{std::back_inserter(y),
    //                                                                         std::back_inserter(y_first_derivative)};

    // ostream iterators
    std::ostream_iterator<double> y_out{std::cout, " ; "};
    std::ostream_iterator<double> y_first_derivative{std::cout, "\n"};

    std::array<std::ostream_iterator<double>, 2> coutput{y_out, y_first_derivative};

    // integration
    runge_kutta::integrate(coutput, inits, differential_equations, 0., 5., 0.1);

    return 0;
}
