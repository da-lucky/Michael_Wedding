#include <algorithm>
#include <array>
#include <functional>
#include <iostream>

namespace runge_kutta
{

template <size_t N, typename T = double>
inline std::array<T, N> _calculate_k(const std::array<std::function<T(T, const std::array<T, N>&)>, N> equations,
                                     const std::array<T, N>& inputs,
                                     T x,
                                     T step)
{
    std::array<T, N> rez{};
    for (size_t i = 0; i < N; ++i)
    {
        rez[i] = equations[i](x, inputs) * step;
    }
    return rez;
}

template <size_t N, typename T = double>
inline std::array<T, N> _prepare_intermediate_input(const std::array<T, N>& init,
                                                    const std::array<T, N>& inputs,
                                                    T multiplier)
{
    std::array<T, N> rez{};

    for (size_t i = 0; i < N; ++i)
    {
        rez[i] = init[i] + (inputs[i] * multiplier);
    }

    return rez;
}

template <size_t N, typename T = double>
inline std::array<T, N> _calculate_outputs(const std::array<T, N>& inputs,
                                           const std::array<T, N>& k1,
                                           const std::array<T, N>& k2,
                                           const std::array<T, N>& k3,
                                           const std::array<T, N>& k4)
{
    std::array<T, N> rez{};

    for (size_t i = 0; i < N; ++i)
    {
        rez[i] = inputs[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }

    return rez;
}

template <typename OutputIterator, size_t N, typename T = double>
inline void _store_output_values(std::array<OutputIterator, N>& out, const std::array<T, N>& vals)
{
    for (size_t i = 0; i < N; ++i)
    {
        out[i] = vals[i];
    }
}

template <typename OutputIterator, size_t N, typename T = double>
void integrate(std::array<OutputIterator, N>& result,
               const std::array<T, N>& initial_values,
               const std::array<std::function<T(T, const std::array<T, N>&)>, N>& equations,
               T interval_start,
               T interval_stop,
               T step)
{
    std::array<T, N> iteration_input_values{initial_values};
    std::array<T, N> iteration_output_values{};
    std::array<T, N> k1_values{};
    std::array<T, N> k2_values{};
    std::array<T, N> k3_values{};
    std::array<T, N> k4_values{};

    _store_output_values(result, initial_values);

    for (T arg{interval_start}; arg <= interval_stop; arg += step)
    {
        k1_values = _calculate_k(equations, iteration_input_values, arg, step);

        k2_values = _calculate_k(
            equations, _prepare_intermediate_input(iteration_input_values, k1_values, 0.5), (arg + step / 2.), step);

        k3_values = _calculate_k(
            equations, _prepare_intermediate_input(iteration_input_values, k2_values, 0.5), (arg + step / 2.), step);

        k4_values = _calculate_k(
            equations, _prepare_intermediate_input(iteration_input_values, k3_values, 1.0), (arg + step), step);

        iteration_output_values =
            _calculate_outputs(iteration_input_values, k1_values, k2_values, k3_values, k4_values);

        iteration_input_values = iteration_output_values;

        _store_output_values(result, iteration_output_values);
    }
}

}  // namespace runge_kutta
