#include <SFML/Graphics.hpp>
#include <iostream>
#include <random>
#include <cmath>

sf::Color from_hsv(float H, float S, float V) {
	float C = S * V;
	float X = C * (1.f - abs(fmod(H / 60.0, 2.) - 1.f));
	float m = V - C;
	float Rs, Gs, Bs;

	if (H >= 0 && H < 60) {
		Rs = C;
		Gs = X;
		Bs = 0;
	}
	else if (H >= 60 && H < 120) {
		Rs = X;
		Gs = C;
		Bs = 0;
	}
	else if (H >= 120 && H < 180) {
		Rs = 0;
		Rs = 0;
		Gs = C;
		Bs = X;
	}
	else if (H >= 180 && H < 240) {
		Rs = 0;
		Gs = X;
		Bs = C;
	}
	else if (H >= 240 && H < 300) {
		Rs = X;
		Gs = 0;
		Bs = C;
	}
	else {
		Rs = C;
		Gs = 0;
		Bs = X;
	}

	sf::Color color;
	color.r = (Rs + m) * 255;
	color.g = (Gs + m) * 255;
	color.b = (Bs + m) * 255;
	return color;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

using fpoint_t = float; // floating point type;

template <class T> T Fmod(T, T) noexcept;
template<> float Fmod<float>(float a, float b) noexcept { return fmodf(a, b);} 
template<> double Fmod<double>(double a, double b) noexcept { return fmod(a, b);} 
template<> long double Fmod<long double>(long double a, long double b) noexcept { return fmodl(a, b);} 











using vec2 = sf::Vector2f;

const fpoint_t pi = 3.14159265358979323846;
const fpoint_t pi2 = 2 * pi;

struct Grid {
    fpoint_t* arr;
    uint size;
    Grid(uint size_) :  arr(new fpoint_t[size_*size_]), size(size_){
        std::random_device r;
        std::default_random_engine e1(r());
        auto distr = std::uniform_real_distribution<fpoint_t>(0, 2 * pi);
        for (uint i = 0u; i < size * size; i++) {
            arr[i] = distr(e1);
        }
    }
    Grid(uint size_, fpoint_t val) :  arr(new fpoint_t[size_*size_]), size(size_){
        for (uint i = 0u; i < size * size; i++) {
            arr[i] = val;
        }
    }


    fpoint_t& operator()(int i, int j) noexcept {
        return arr[i * size + j];
    }
    ~Grid() {
        delete[] arr;
    }
};

struct MagnetSystem_T2 {
    Grid angles;
    Grid vels;
    uint size;
    fpoint_t alpha;
    fpoint_t damp;
    double time;
    MagnetSystem_T2(uint size_, fpoint_t alpha_, fpoint_t damp_) : angles(size_), vels(size_, 0), size(size_), alpha(alpha_), damp(damp_), time(0) {}

    fpoint_t& operator()(uint i, uint j) noexcept {
        return angles(i, j);
    }
    void Step(double dt) noexcept{
        time += dt;
        for (uint i = 0u; i < size; i++) {
            for (uint j = 0u; j < size; j++) {
                fpoint_t& phi = angles(i, j);
                fpoint_t& vel = vels(i, j);

                phi = phi + vel * dt;  

                fpoint_t accel = 0.f;
                accel += sin(phi - Angle_0_dec(i, j));
                accel += sin(phi - Angle_0_incr(i, j));
                accel += sin(phi - Angle_dec_0(i, j));
                accel += sin(phi - Angle_incr_0(i, j));
                accel *= alpha;

                vel += accel * dt;
                vel -=  damp * vel;
            }
        }
    }

    // Two-dimensional iterators
    fpoint_t Angle_dec_0(uint i, uint j) noexcept {
        if  (i != 0u) { return angles(i-1u, j); }
        return angles(size - 1, j);
    } 
    fpoint_t Angle_incr_0(uint i, uint j) noexcept {
        if  (i != size - 1u) {  return angles(i+1u, j); }
        return angles(0, j);
    } 
    fpoint_t Angle_0_dec(uint i, uint j) noexcept {
        if  (j != 0u) { return angles(i, j-1u); }
        return angles(i, size - 1);
    } 
    fpoint_t Angle_0_incr(uint i, uint j) noexcept {
        if  (j != size - 1u) { return angles(i, j+1u); }
        return angles(i, 0);
    }
    // Avoid overflow, increase precision
    void NormalizeAngles() noexcept {
        for (uint i = 0u; i < size; i++) {
            for (uint j = 0u; j < size; j++) {
                fpoint_t& phi = angles(i, j);
                phi = Fmod(phi, pi);  
            }
        }
    }
    // Do many steps at one
    void BatchStep(uint count, double dt) noexcept {
        for (uint i = 0; i < count; i++) {
            Step(dt);
        }
        NormalizeAngles();
    };
    template <uint count>
    void InlineBatchStep(double dt) noexcept {
        for (uint i = 0; i < count; i++) {
            Step(dt);
        }
        NormalizeAngles();
    }
    long double EvalEnergy() noexcept {
        long double energy = 0;
        for (uint i = 0u; i < size; i++) {
            for (uint j = 0u; j < size; j++) {
                fpoint_t& phi = angles(i, j);
                fpoint_t& vel = vels(i, j);

                energy += cos(phi - Angle_0_dec(i, j)) * alpha / 2;
                energy += cos(phi - Angle_0_incr(i, j)) * alpha / 2;
                energy += cos(phi - Angle_dec_0(i, j)) * alpha / 2;
                energy += cos(phi - Angle_incr_0(i, j)) * alpha / 2;
                energy += vel * vel / 2.;
            }
        }
        return energy * alpha;
    }

};

// Drawing mode
enum What_Show { 
    show_none, // black screen
    show_colors, // angles drawn as hues
    show_smooth_colors, // monochrome mode
    show_heat, // show velocity aka heat, log
    COUNT
};

int main()
{
    const int step = 5; // size of a pixel
    const int pixsize = 800; // screen size
    What_Show what_show = show_colors;
    // create the window
    sf::RenderWindow window(sf::VideoMode(pixsize, pixsize), "Simulation");

    sf::RectangleShape rectangle(sf::Vector2f(step, step));
    MagnetSystem_T2 grid(pixsize / step, -1., 0.005);


    // run the program as long as the window is open
    while (window.isOpen())
    {
        // end the previous frame
        window.display();

        // check all the window's events that were triggered since the last iteration of the loop
        sf::Event event;
        while (window.pollEvent(event))
        {
            switch(event.type) {
                case(sf::Event::Closed) : {
                    window.close();
                    break;
                }
                case(sf::Event::KeyPressed) :  {
                    switch(event.key.code) {
                        case(sf::Keyboard::N) : {
                            what_show = show_none;
                            break;
                        }
                        case(sf::Keyboard::C) : {
                            what_show = show_colors;
                            break;
                        }
                        case(sf::Keyboard::H) : {
                            what_show = show_heat;
                            break;
                        }
                        case(sf::Keyboard::V) : {
                            what_show = show_smooth_colors;
                            break;
                        }
                        case(sf::Keyboard::E) : {
                            std::cout << grid.EvalEnergy() << std::endl;
                            break;
                        }
                        case(sf::Keyboard::D) : {
                            std::cin >> grid.damp;
                            break;
                        }
                        default : {
                            break;
                        }
                    }
                    break;
                }
                case (sf::Event::MouseButtonPressed) : {
                    if (event.mouseButton.button == sf::Mouse::Left) {
                        uint i = event.mouseButton.x / step;
                        uint j = event.mouseButton.y / step;
                        if (i < grid.size && j < grid.size) {
                            for (int k = -10; k < 10; k++) {
                                for (int l = -10; l < 10; l++) {
                                    grid.vels(i + k, j + l)  += 0.5;
                                }
                            }
                        }
                    } else if (event.mouseButton.button == sf::Mouse::Right) {
                        uint i = event.mouseButton.x / step;
                        uint j = event.mouseButton.y / step;
                        if (i < grid.size && j < grid.size) {
                            for (int k = -10; k < 10; k++) {
                                for (int l = -10; l < 10; l++) {
                                    grid.angles(i + k, j + l) += pi / 2.;
                                }
                            }
                        }
                    }
                    break;
                }
                default : {
                    break;
                }
            }
        }


        
        window.clear(sf::Color::Black);
        grid.InlineBatchStep<100>(0.0005);

        if (what_show == show_none) {
            continue;
        } else if (what_show == show_colors) {
            for (uint i = 0u; i < grid.size; i++) {
                for (uint j = 0u; j < grid.size; j++) {
                    float hue = 180.f / pi * fmod(grid(i, j) + pi2, pi2);
                    rectangle.setFillColor(from_hsv(hue, 1, 0.8));
                    rectangle.setPosition(sf::Vector2f(float(i * step), float(j * step)));
                    window.draw(rectangle);
                }
            }
        } else if (what_show == show_smooth_colors) {
            for (uint i = 0u; i < grid.size; i++) {
                for (uint j = 0u; j < grid.size; j++) {
                    float hue = sin(grid(i, j));
                    rectangle.setFillColor(sf::Color(127.f + 126.f * hue, 0, 0));
                    rectangle.setPosition(sf::Vector2f(float(i * step), float(j * step)));
                    window.draw(rectangle);
                }
            }

        } else if  (what_show == show_heat) {
            for (uint i = 0u; i < grid.size; i++) {
                for (uint j = 0u; j < grid.size; j++) {
                    float param = 0.1 * std::max(log(grid.vels(i, j) * grid.vels(i, j)) / 2 + 10, 0.);
                    rectangle.setFillColor(sf::Color(255, 0, 0, 255 * param));
                    rectangle.setPosition(sf::Vector2f(float(i * step), float(j * step)));
                    window.draw(rectangle);
                }
            }
        }
    }

    return 0;
}