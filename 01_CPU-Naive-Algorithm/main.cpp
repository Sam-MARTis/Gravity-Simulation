#include <iostream>
#include <stdlib.h>
#include <random>
#include <SFML/Graphics.hpp>


// Simulation properties
#define PARTICLE_COUNT 100
#define RAD 2.0f
#define dRad 1.0f


// Physics properties
#define G 1000.0f
#define RAD_MASS_CONSTANT 100.0f

// //Memory objects

struct Particle{
    float x;
    float y;
    float rad;
    float mass;
};



//Helper functions
float randf(float lb = 0.0f, float ub = 1.0f){
    return lb + (ub-lb)*(((float)rand())/((float)RAND_MAX));
}




// Supplementary functions
void initialize_particles(float* pos, float * mas, sf::CircleShape* shapes){

}




int main()
{
    return 0;
}