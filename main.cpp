#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <bitset>

using interval_t = std::pair <float, float>;
using arr_t = std::vector <float>;
using arr_it = arr_t::iterator;

float calculate_polunom(float x, const arr_t& a)
{
	auto ans = 0.0f;
	for (auto i = 0U; i < std::size(a); i++)
	{
		ans += std::powf(x, i) * a[i];
	}
	return ans;
}

template < unsigned N >
std::vector < unsigned short > find_number_of_sign_change(float x_0, float x_last, const std::vector <arr_t>& sturm_series)
{
	std::vector < unsigned short > number_of_sign_change(N, 0Ull);
	std::bitset< N > prev_sign;

	float y;
	auto step = (x_last - x_0) / N;
	auto x = x_0;
	for (auto i = 0U; i < N; i++)
	{
		y = calculate_polunom(x, sturm_series[0]);
		if (y >= 0.0f)
		{
			prev_sign[i] = true;
		}
		else
		{
			prev_sign[i] = false;
		}
		x += step;
	}
	for (auto n = 1U; n < std::size(sturm_series); n++)
	{
		x = x_0;
		for (auto i = 0U; i < N; i++)
		{
			y = calculate_polunom(x, sturm_series[n]);
			if ((y > 0.0f && !prev_sign[i]) || (y < 0.0f && prev_sign[i]))
			{
				number_of_sign_change[i]++;

				if (y >= 0.0f)
				{
					prev_sign[i] = true;
				}
				else
				{
					prev_sign[i] = false;
				}
			}
			x += step;
		}
	}
	return number_of_sign_change;
}

arr_t polynomial_remainder(const arr_t& a, const arr_t& b)
{
	auto m = std::size(b) - 1;
	auto n = std::size(a) - 1;

	if (m > n)
	{
		throw std::logic_error("n should be greater than m");
	}
	if (b[m] == 0.0f)
	{
		throw std::logic_error("b_m should be greater than 0");
	}

	arr_t c(n + 1);
	std::copy(std::begin(a), std::end(a), std::begin(c));

	float k_i;
	for (auto i = 0U; i < n - m + 1U; i++)
	{
		k_i = c[n - i] / b[m];
		for (auto j = i; j <= m + i; j++)
		{
			c[n - j] -= k_i * b[m - j + i];
		}
		c.pop_back();
	} 
	while (c[std::size(c) - 1] == 0.0)
	{
		c.pop_back();
	}
 	return c;
}

arr_t polynomial_deritivate(const arr_t& a)
{
	// a[n]Z^n + a[n - 1]Z^(n-1) + ... + a[1]Z + a[0]
	// b[n-1] = a[n]*n
	// ...
	// b[1] = a[2] * 2
	// b[0] = a[1] * 1

	arr_t ans(std::size(a) - 1U);

	for (auto i = 1U; i < std::size(a); i++)
	{
		ans[i - 1U] = a[i] * i;
	}
	return ans;
}

void build_sturm_series(const arr_t& a, std::vector <arr_t>& output)
{
	output.emplace_back(a);

	arr_t der_a = polynomial_deritivate(a);
	output.emplace_back(der_a);

	auto i = 3U;
	arr_t rem = polynomial_remainder(output[0], output[1]);
	while (std::size(rem) > 1)
	{
		std::transform(std::begin(rem), std::end(rem), std::begin(rem), [](auto x) {return -x; });
		output.emplace_back(rem);
		rem = polynomial_remainder(output[i - 2U], output[i - 1U]);
		i++;
	}
	std::transform(std::begin(rem), std::end(rem), std::begin(rem), [](auto x) {return -x; });
	output.emplace_back(rem);
}

void localize_sturm(arr_t& a, std::vector <interval_t> & localized_intervals, arr_t & roots)
{
	arr_t abs_a(std::size(a));
	std::transform(std::begin(a), std::end(a), std::begin(abs_a), [](auto x) { return std::abs(x); });
	auto A = *std::max_element(std::begin(abs_a), std::end(abs_a) - 1);
	auto B = *std::max_element(std::begin(abs_a) + 1, std::end(abs_a));
	interval_t init_positive_interval = { *std::begin(abs_a) / (*std::begin(abs_a) + B) , 1.0f + A / *(std::end(abs_a) - 1U)};

	std::vector < arr_t > sturm_series;
	build_sturm_series(a, sturm_series);
	
	const auto T_n = 100ULL; 
	const auto x_0 = - init_positive_interval.second;
	const auto x_last = init_positive_interval.second;
	const auto step = (x_last - x_0) / T_n;
	auto number_of_sign_change = find_number_of_sign_change< T_n >(x_0, x_last, sturm_series);
	
	for (auto i = 1U; i < std::size(number_of_sign_change); i++)
	{
		if (number_of_sign_change[i - 1U] - number_of_sign_change[i] >= 1U)
		{
			localized_intervals.push_back(std::make_pair(x_0 + step * (i - 1U), x_0 + step * i));
		}
	}
}

float accure_root(const arr_t & a, interval_t interval, const float eps = 1.0e-6f)
{
	auto l = interval.first;
	auto r = interval.second;
	auto c = r + (l - r) * 0.5f;
	while (r - l > eps)
	{
		if (calculate_polunom(l, a) * calculate_polunom(c, a) < 0)
		{
			r = c;
		}
		else
		{
			l = c;
		}
		c = r + (l - r) * 0.5f;
	}
	return c;
}

int main(int argc, char** argv)
{
	// initialize parameters
	constexpr float gamma0 = 5.f / 3.f;
	constexpr float dens0 = 1.0e-5f;
	constexpr float U0 = 0.0f;
	constexpr float P0 = 3.848e3f;
	constexpr float gamma3 = 5.f / 3.f;
	constexpr float C3 = 2.53248e4f;
	constexpr float U3 = 0.0f;
	constexpr float P3 = 3.04e9f;

	// compute 
	constexpr float dens3 = gamma3 * P3 / C3 / C3;
	constexpr float alpha0 = (gamma0 + 1.0f) / (gamma0 - 1.0f);
	constexpr std::size_t n = 2.0f * gamma3 / (gamma3 - 1.0f); // should be integer
	const float mu = (U3 - U0) * std::sqrtf((gamma0 - 1.0f) * dens0 / 2.0f / P0);
	const float nu = 2.0f / (gamma3 - 1.0f) * std::sqrtf(gamma3 * (gamma0 - 1.0f) * P3 * dens0 / 2.0f / P0 / dens3);
	constexpr float X = P3 / P0;

	// a[n]Z^n + a[n - 1]Z^(n-1) + ... + a[1]Z + a[0] 
	arr_t a(2U * n + 1U, 0.0f);
	a[2ULL * n] = X * X;
	a[n + 2ULL] = - alpha0 * nu * nu * X;
	a[n + 1ULL] = 2.0f * alpha0 * nu * (nu + mu) * X;
	a[n] = - (2.0f + (nu + mu) * (nu + mu) * alpha0) * X;
	a[2ULL] = - nu * nu;
	a[1ULL] = 2.0f * nu * (nu + mu);
	a[0ULL] = 1.0f - (nu + mu) * (nu + mu);

	arr_t roots;

	//std::cout << calculate_polunom(4.0f, b) << '\n';
	std::vector <interval_t> localized_intervals;
	localize_sturm(a, localized_intervals, roots);

	for (const auto& interval : localized_intervals)
	{
		if (calculate_polunom(interval.first, a) * calculate_polunom(interval.second, a) < 0.0f)
		{
			auto root = accure_root(a, interval);
			std::cout << '(' << interval.first << ' ' << interval.second << ")		--> " << root << '\n';
			roots.push_back(root);
		}
	}

	auto Z = roots[std::size(roots) - 1ULL];
	auto P1 = P3 * std::powf(Z, n);
	auto dens1 = dens0 * ((gamma0 - 1.0f) * P0 + (gamma0 + 1.0f) * P1) / ((gamma0 + 1.0f) * P0 + (gamma0 - 1.0f) * P1);
	auto U1_U0 = std::sqrt( (dens1 + dens0) * (P1 - P0) / dens1 / dens0);
	auto D0 = (P1 - P0) / dens0 / U1_U0 - U0;
	std::cout << "D_0 = " << D0 << std::endl;
}