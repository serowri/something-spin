#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <limits>
#include <cstdint>

static constexpr double pi = 3.141592653589793;
static constexpr double a = -pi/16, b = pi/32, c = pi/32; // angle parameters (x, y ,z) = (a, b, c).
static constexpr double p_z = 3.5, f_z = 6.0; // projection surface z, focus point z. (x, y) = (0, 0)
static constexpr double N_min = -5.0; // instead of numerical_limits

static constexpr double n_pixels = 61.0, width = 2.5; // terminal output string numbers, projection surface width
static constexpr double ratio = n_pixels / width;
static constexpr uint64_t surface_pnums = 40, xyz = 3; // surface pixel numbers
static constexpr double side_len = 2.0; // side lenght of cude
static constexpr uint64_t n_millis = 100; // update span (ms)

// Note: font size (row, column): (1.2, 0.6)(default) -> (0.8, 0.8).

std::vector<std::vector<std::pair<std::string, double>>> p_surface(n_pixels, std::vector<std::pair<std::string, double>>(n_pixels));

class Surface {
	std::vector<std::vector<std::vector<double>>> s{xyz, std::vector<std::vector<double>>(surface_pnums, std::vector<double>(surface_pnums))};
	std::string surface_code;
	double dx, dy, dz;
	double b_dx, b_dy, b_dz;

public:
	Surface (std::vector<double> ltop, std::vector<std::pair<bool, double>> diff, std::string code): 
		b_dx(diff.at(0).first), b_dy(diff.at(1).first), b_dz(diff.at(2).first), dx(diff.at(0).second), dy(diff.at(1).second), dz(diff.at(2).second), surface_code(code) {
		
		// s init
		double x = ltop.at(0), y = ltop.at(1), z = ltop.at(2);
		for (uint64_t i = 0; i < surface_pnums; ++i) {
			for (uint64_t j = 0; j < surface_pnums; ++j) {
				s.at(0).at(i).at(j) = x;
				s.at(1).at(i).at(j) = y;
				s.at(2).at(i).at(j) = z;
				if (!b_dx) {
					z += dz;
				}else if (!b_dy) {
					x += dx;
				}else if (!b_dz) {
					y += dy;
				}

			}

			if (!b_dx) {
				y += dy;
				z = ltop.at(2);
			}else if (!b_dy) {
				z += dz;
				x = ltop.at(0);
			}else if (!b_dz) {
				x += dx;
				y = ltop.at(1);
			}
		}
		// s init
		
	}

	void cal_xyz () { 
		// calculate next step coord :(new_x, new_y, new_z)^T = Rz(c)*Ry(b)*Rx(a)*(x, y ,z)^T
		
		for (uint64_t i = 0; i < surface_pnums; ++i) {
			for (uint64_t j = 0; j < surface_pnums; ++j) {
				double x = s.at(0).at(i).at(j); 
			    double y = s.at(1).at(i).at(j); 
				double z = s.at(2).at(i).at(j);
				double new_x = x; 
				double new_y = (std::cos(a)*y - std::sin(a)*z);
				double new_z = (std::sin(a)*y + std::cos(a)*z);
				x = new_x;
				y = new_y;
				z = new_z;
				new_x = (std::cos(b)*x + std::sin(b)*z);
				new_y = y;
				new_z = (-std::sin(b)*x + std::cos(b)*z);
				x = new_x;
				y = new_y;
				z = new_z;
				new_x = (std::cos(c)*x - std::sin(c)*y);
				new_y = (std::sin(c)*x + std::cos(c)*y);
				new_z = z;
				s.at(0).at(i).at(j) = new_x;
				s.at(1).at(i).at(j) = new_y; 
				s.at(2).at(i).at(j) = new_z;
			}
		}
		return;
	}

	std::tuple<double, double, double> cal_ps_xyz (double x, double y, double z) {
		// calculate projection surface coord
		double p_x = (p_z - z) * (x / (f_z - z));
		double p_y = (p_z - z) * (y / (f_z - z));
		return {p_x, p_y, z}; // (p_x, p_y) & z
	}

	void update_surface() {
		// translate projection coord to vector index and update its surface code
		double ps_x, ps_y, ps_z;
		double x_start = -(n_pixels / 2) / ratio;
		double y_start = -x_start;
		for (uint64_t i = 0; i < surface_pnums; ++i) {
			for (uint64_t j = 0; j < surface_pnums; ++j) {
				double x = s.at(0).at(i).at(j);
				double y = s.at(1).at(i).at(j); 
				double z = s.at(2).at(i).at(j);
				std::tie(ps_x, ps_y, ps_z) = this->cal_ps_xyz(x, y, z);
				uint64_t target_i = static_cast<uint64_t>(std::floor(ratio * (y_start - ps_y))); // y~i floor?
				uint64_t target_j = static_cast<uint64_t>(std::floor(ratio * (ps_x - x_start))); // x~j
				if (target_i < 0 || target_i > n_pixels-1) {
					std::cout << "out of range: (i, j)"  << "(" << target_i << ", " << target_j << ")" << std::endl;
					return;
				}
				if (target_j < 0 || target_j > n_pixels-1) {
					std::cout << "out of range: (i, j)"  << "(" << target_i << ", " << target_j << ")" << std::endl;
					return;
				}
				if (z > p_surface.at(target_i).at(target_j).second) {
					p_surface.at(target_i).at(target_j).first = surface_code;
					p_surface.at(target_i).at(target_j).second = z;
				}
			}
		}
		return;
	}
};



// all
void reset_surface () {
	// reset surface code and depth( = z)
	for (uint64_t i = 0; i < p_surface.size(); ++i) {
		for (uint64_t j = 0; j < p_surface.at(0).size(); ++j) {
			p_surface.at(i).at(j).first = " ";
			p_surface.at(i).at(j).second = N_min;
		}
	}
	return;
}

// all
void print_cube () { 
	// clear and print surface code
	std::cout << "\033[2J\033[1;1H";
	for (uint64_t i = 0; i < p_surface.size(); ++i) {
		for (uint64_t j = 0; j < p_surface.at(0).size(); ++j) {
			std::cout << p_surface.at(i).at(j).first;
		}
		std::cout << std::endl;
	}
	return;
}


int main () {
	
	std::cout << "\033[2J\033[1;1H";
	// init 6 surfaces
	double d = side_len / static_cast<double>(surface_pnums); 
	Surface s_A({1.0, 1.0, 1.0}, {{false, 0.0}, {true, -d}, {true, -d}}, "|");
	Surface s_B({-1.0, 1.0, 1.0}, {{true, d}, {true, -d}, {false, 0.0}}, "#");
	Surface s_C({-1.0, 1.0, -1.0}, {{false, 0.0}, {true, -d}, {true, d}}, "_");
	Surface s_D({1.0, 1.0, -1.0}, {{true, -d}, {true, -d}, {false, 0.0}}, "/");
	Surface s_E({-1.0, 1.0, 1.0}, {{true, d}, {false, 0.0}, {false, -d}}, "+");
	Surface s_F({-1.0, -1.0, -1.0}, {{true, d}, {false, 0.0}, {true, d}}, "-");
	// init 6 surfaces

	while (1) {
		std::this_thread::sleep_for(std::chrono::milliseconds(n_millis));

		s_A.cal_xyz();
		s_B.cal_xyz();
		s_C.cal_xyz();
		s_D.cal_xyz();
		s_E.cal_xyz();
		s_F.cal_xyz();
		
		reset_surface();

		s_A.update_surface();
		s_B.update_surface();
		s_C.update_surface();
		s_D.update_surface();
		s_E.update_surface();
		s_F.update_surface();
		
		print_cube();
	}
	return 0;
}

