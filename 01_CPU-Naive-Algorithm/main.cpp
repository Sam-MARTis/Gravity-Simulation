#include <iostream>
#include <stdlib.h>
#include <random>
#include <SFML/Graphics.hpp>


// Simulation properties
#define PARTICLE_COUNT 100
#define RAD 2.0f
#define dRad 0.0f
#define XMin 000.0f
#define XMax 800.0f
#define YMin 000.0f
#define YMax 600.0f
#define MAGNIFICATION 0.001f
#define MASS 10.0f

#define SUBSTEPS 501



#define SUBSTEPS_DT 0.01f


// Physics properties
#define G 1e-7f
#define RAD_MASS_CONSTANT 10.0f
#define MAX_FORCE 1e-3f

// //Memory objects

struct Particle{
    float x;
    float y;
    float xprev;
    float yprev;
    float vx;
    float vy;
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
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        sf::CircleShape& shape = shapes[i];
        particle.x = randf(XMin*MAGNIFICATION, XMax*MAGNIFICATION);
        particle.y = randf(YMin*MAGNIFICATION, YMax*MAGNIFICATION);
        particle.xprev = particle.x;
        particle.yprev = particle.y;
        particle.vx = 0.0f;
        particle.vy = 0.0f;
        particle.rad = RAD + randf(-dRad, dRad);
    
        // particle.mass = RAD_MASS_CONSTANT * particle.rad * particle.rad;
        particle.mass = MASS;
        shape.setRadius(particle.rad);
        shape.setOrigin(particle.rad, particle.rad);
        shape.setPosition(particle.x, particle.y);
        shape.setFillColor(sf::Color::White);
    }
    particles[count-1].x = MAGNIFICATION* (XMax + XMin) / 2.0f;
    particles[count-1].y = MAGNIFICATION* (YMax + YMin) / 2.0f;
    particles[count-1].xprev = particles[count-1].x;
    particles[count-1].yprev = particles[count-1].y;
    particles[count-1].rad = 3*RAD;
    particles[count-1].mass = 1000.0f;
    particles[count-1].vx = 0.0f;
    particles[count-1].vy = 0.0f;
    shapes[count-1].setRadius(particles[count-1].rad);
    shapes[count-1].setOrigin(particles[count-1].rad, particles[count-1].rad);
    shapes[count-1].setPosition(particles[count-1].x, particles[count-1].y);
    shapes[count-1].setFillColor(sf::Color::Red);


    std::cout << "Particles initialized." << std::endl;
}


void step_gravity(Particle* particles, int count){
    // Verlet integration for gravity between particles
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        float fx = 0.0f;
        float fy = 0.0f;

        for(int j=0; j<count; j++){
            if(i == j) continue; // Skip self

            Particle& other = particles[j];
            float dx = (other.x - particle.x)*100;
            float dy = (other.y - particle.y)*100;
            float dist_sq = dx * dx + dy * dy;


            if(dist_sq < 1e-6f) continue; // Avoid division by zero

            const float dist = sqrtf(dist_sq);
            float force = G * (particle.mass * other.mass) / dist_sq;
            force = force > MAX_FORCE ? MAX_FORCE : force; // Limit force to MAX_FORCE
            fx += force * (dx / dist);
            fy += force * (dy / dist);

        }
        // if(fx*fx + fy*fy > MAX_FORCE * MAX_FORCE){
        //     float norm = sqrtf(fx*fx + fy*fy);
        //     fx = (fx / norm) * MAX_FORCE;
        //     fy = (fy / norm) * MAX_FORCE;
        // }
        // Verlet integration
        const float dt = SUBSTEPS_DT;
        // float xnew = particle.x + (particle.x - particle.xprev) + fx * dt * dt / particle.mass;
        // float ynew = particle.y + (particle.y - particle.yprev) + fy * dt * dt / particle.mass;
        // particle.xprev = particle.x;
        // particle.yprev = particle.y;
        // particle.x = xnew;
        // particle.y = ynew;
        particle.vx += fx * dt / particle.mass;
        particle.vy += fy * dt / particle.mass;
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;


    }

}

void sync_shapes(Particle* particles, sf::CircleShape* shapes, int count){
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


    sf::RenderWindow window(sf::VideoMode(800, 600), "Particle Simulation");
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