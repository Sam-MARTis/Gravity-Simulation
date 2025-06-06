#include <iostream>
#include <stdlib.h>
#include <random>
#include <SFML/Graphics.hpp>
#include <omp.h>

/*
g++ main.cpp -o grav -lsfml-graphics -lsfml-window -lsfml-system -fopenmp
*/
// Simulation properties
#define PARTICLE_COUNT 3000
#define RAD 1.0f
#define dRad 0.0f
#define XMin 000.0f
#define XMax 800.0f
#define YMin 000.0f
#define YMax 800.0f
#define MAGNIFICATION 1.0f
#define MASS 10.0f
#define CENTER_MASS 10000.0f

#define SUBSTEPS 10



#define SUBSTEPS_DT 0.01f


// Physics properties
#define G 0.10f
#define RAD_MASS_CONSTANT 10.0f
#define MAX_FORCE 1000.0f

// //Memory objects

struct Particle{
    float x;
    float y;
    float vx;
    float vy;
    float ax;
    float ay;
    float rad;
    float mass;
};

Particle *particles = new Particle[PARTICLE_COUNT];
sf::CircleShape *pshape = new sf::CircleShape[PARTICLE_COUNT];


//Helper functions
float randf(float lb = 0.0f, float ub = 1.0f){
    return lb + (ub-lb)*(((float)rand())/((float)RAND_MAX));
}




// Supplementary functions
void initialize_particles(Particle* particles, sf::CircleShape* shapes, int count){
    const float centerX = (XMax + XMin) / 2.0f;
    const float centerY = (YMax + YMin) / 2.0f;
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        sf::CircleShape& shape = shapes[i];
        // particle.x = randf(XMin*MAGNIFICATION, XMax*MAGNIFICATION);
        // particle.y = randf(YMin*MAGNIFICATION, YMax*MAGNIFICATION);
        // Initialize velocity to be in stable circular motion. centripital == gravitational force
        float angle = randf(0.0f, 2.0f * 3.14159265358979323846f);
        const float radius = MAGNIFICATION * randf(0.0f, 0.7f* (XMax - XMin) / 2.0f); // Use half the width as radius
        particle.x = MAGNIFICATION * ((XMax + XMin) / 2.0f + radius * cosf(angle));
        particle.y = MAGNIFICATION * ((YMax + YMin) / 2.0f + radius * sinf(angle));
        // particle.x = MAGNIFICATION * (XMax + XMin) / 2.0f + RAD * cosf(angle);
        // particle.y = MAGNIFICATION * (YMax + YMin) / 2.0f + RAD * sinf(angle);
        float speed =  sqrtf(G * CENTER_MASS / (radius)); // Speed for stable circular motion
        particle.vx = speed * sinf(angle);
        particle.vy = -speed * cosf(angle);
        
        particle.ax = 0.0f;
        particle.ay = 0.0f;
        particle.rad = RAD + randf(-dRad, dRad);
    
        // particle.mass = RAD_MASS_CONSTANT * particle.rad * particle.rad;
        particle.mass = MASS;
        shape.setRadius(particle.rad);
        shape.setOrigin(particle.rad, particle.rad);
        shape.setPosition(particle.x, particle.y);
        shape.setFillColor(sf::Color::White);
    }
    particles[0].x = MAGNIFICATION* (XMax + XMin) / 2.0f;
    particles[0].y = MAGNIFICATION* (YMax + YMin) / 2.0f;
    particles[0].rad = 3*RAD;
    particles[0].mass = CENTER_MASS;
    particles[0].vx = 0.0f;
    particles[0].vy = 0.0f;
    particles[0].ax = 0.0f;
    particles[0].ay = 0.0f;
    shapes[0].setRadius(particles[0].rad);
    shapes[0].setOrigin(particles[0].rad, particles[0].rad);
    shapes[0].setPosition(particles[0].x, particles[0].y);
    shapes[0].setFillColor(sf::Color::Red);


    std::cout << "Particles initialized." << std::endl;
}


void step_gravity(Particle* particles, int count){
    // Verlet integration for gravity between particles
    #pragma omp parallel for schedule(static)
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        float fx = 0.0f;
        float fy = 0.0f;
        #pragma omp parralel for
        for(int j=0; j<count; j++){
            if(i == j) continue; // Skip self

            Particle& other = particles[j];
            float dx = (other.x - particle.x);
            float dy = (other.y - particle.y);
            float dist_sq = dx * dx + dy * dy;


            if(dist_sq < 1e-2f) continue; // Avoid division by zero

            const float dist = sqrtf(dist_sq);
            float force = G * (particle.mass * other.mass) / dist_sq;
            force = force > MAX_FORCE ? 0 : force; // Limit force to MAX_FORCE
            particle.ax += force * (dx / dist)/particle.mass;
            particle.ay += force * (dy / dist)/particle.mass;
            // other.ax -= force * (dx / dist)/other.mass;
            // other.ay -= force * (dy / dist)/other.mass;

        }

    }
    // Update positions and velocities
    #pragma omp parallel for
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        particle.x += particle.vx * SUBSTEPS_DT*0.5f;
        particle.y += particle.vy * SUBSTEPS_DT*0.5f;
        // if(i==100) std::cout << "Particle " << i << " position: (" << particle.x/MAGNIFICATION << ", " << particle.y/MAGNIFICATION << ")" << std::endl;
        particle.vx += particle.ax * SUBSTEPS_DT;
        particle.vy += particle.ay * SUBSTEPS_DT;
        particle.x += particle.vx * SUBSTEPS_DT*0.5f;
        particle.y += particle.vy * SUBSTEPS_DT*0.5f;
        particle.ax = 0.0f; 
        particle.ay = 0.0f;
    }

}

void sync_shapes(Particle* particles, sf::CircleShape* shapes, int count){
    #pragma omp parallel for
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        sf::CircleShape& shape = shapes[i];
        shape.setPosition(particle.x/MAGNIFICATION, particle.y/MAGNIFICATION);
    }
}


int main()
{


    // Initialize random seed
    srand(42);
    initialize_particles(particles, pshape, PARTICLE_COUNT);


    sf::RenderWindow window(sf::VideoMode(XMax-XMin, YMax-YMin), "Particle Simulation");
    window.setFramerateLimit(60);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::Black);
        
        for (int i =0; i< SUBSTEPS; i++){
            // Perform substeps
            step_gravity(particles, PARTICLE_COUNT);
        }
        sync_shapes(particles, pshape, PARTICLE_COUNT);

        for(int i=0; i<PARTICLE_COUNT; i++){
            window.draw(pshape[i]);
        }

        window.display();
    }
    return 0;
}