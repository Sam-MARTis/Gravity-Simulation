#include <iostream>
#include <stdlib.h>
#include <random>
#include <SFML/Graphics.hpp>


// Simulation properties
#define PARTICLE_COUNT 1000
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




#define DEBUG



int current_available_index = PARTICLE_COUNT;

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

struct Node{
    
    bool valid;
    int leaves[4];
    float centerX;
    float centerY;
    float sideLength;
    
    Node() : valid(false), leaves{-1, -1, -1, -1}, centerX(0.0f), centerY(0.0f), sideLength(0.0f) {}
};



Particle *particles = new Particle[PARTICLE_COUNT];
sf::CircleShape *pshape = new sf::CircleShape[PARTICLE_COUNT];


//Helper functions
float randf(float lb = 0.0f, float ub = 1.0f){
    return lb + (ub-lb)*(((float)rand())/((float)RAND_MAX));
}
void print(auto message){
    std::cout << message << std::endl;
}




// Supplementary functions


void init_particles(Particle* particles, sf::CircleShape* shapes, int count){
    const float centerX = (XMax + XMin) / 2.0f;
    const float centerY = (YMax + YMin) / 2.0f;
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        sf::CircleShape& shape = shapes[i];
        // particle.x = randf(XMin*MAGNIFICATION, XMax*MAGNIFICATION);
        // particle.y = randf(YMin*MAGNIFICATION, YMax*MAGNIFICATION);
        // Initialize velocity to be in stable circular motion. centripital == gravitational force
        float angle = randf(0.0f, 2.0f * 3.14159265358979323846f);
        const float radius = MAGNIFICATION * randf(0.0f, 0.5f* (XMax - XMin) / 2.0f); // Use half the width as radius
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
void compute_forces_naive(Particle* particles, int count){
    // Compute forces between particles
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        float fx = 0.0f;
        float fy = 0.0f;

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

        }

    }
}

// void insert_point(int *tree_array, int)


void get_bounding_box(Particle* particles, int count, float* bounds){
    float minX = particles[0].x;
    float minY = particles[0].y;
    float maxX = particles[0].x;
    float maxY = particles[0].y;
    
    for(int i=1; i<count; i++){
        Particle& particle = particles[i];
        if(particle.x < minX) minX = particle.x;
        if(particle.y < minY) minY = particle.y;
        if(particle.x > maxX) maxX = particle.x;
        if(particle.y > maxY) maxY = particle.y;
    }
    bounds[0] =( minX + maxX) * 0.5f;
    bounds[1] = (minY + maxY) * 0.5f;
    bounds[2] = ((maxX - minX)> (maxY - minY)? (maxX - minX) : (maxY - minY))*0.5f; //Square grid

    // bounds[0] = minX;
    // bounds[1] = minY;
    // bounds[2] = maxX;
    // bounds[3] = maxY;
}


void insert_particle(Particle* particles, const int pidx, Node* tree_array, const int index){
    Node& node = tree_array[index];
    #ifdef DEBUG
    if(!node.valid) print("Node is not valid, cannot insert particle.");
    #endif
    Particle& particle = particles[pidx];
    int insert_index = (particle.x > node.centerX) + 2 * (particle.y > node.centerY);
    if(node.leaves[insert_index] == -1){
        node.leaves[insert_index] = pidx;
        return;
    }else if (node.leaves[insert_index] >= PARTICLE_COUNT) //It has a node. Go to the node pointed to
    {
        insert_particle(particles, pidx, tree_array, node.leaves[insert_index] - PARTICLE_COUNT);
        return;
    }

    else if (node.leaves[insert_index]< PARTICLE_COUNT)
    {
        // Node is a leaf, we need to split it
        int old_particle_idx = node.leaves[insert_index];
        int new_node_idx = current_available_index++;
        node.leaves[insert_index] = new_node_idx + PARTICLE_COUNT ; // Update the current node to point to the new node
        Node& new_node = tree_array[new_node_idx];
        new_node.valid = true;
        // if(insert_index == 0){
        //     new_node.centerY = node.centerY - node.sideLength * 0.5f; // Upper half
        //     new_node.centerX = node.centerX - node.sideLength * 0.5f; // Right half
        // }
        // if(insert_index == 1){
        //     new_node.centerY = node.centerY - node.sideLength * 0.5f; // Lower half
        //     new_node.centerX = node.centerX + node.sideLength * 0.5f; // Left half
        // }
        // if(insert_index == 2){
        //     new_node.centerY = node.centerY + node.sideLength * 0.5f; // Upper half
        //     new_node.centerX = node.centerX - node.sideLength * 0.5f; // Right half
        // }
        // if(insert_index == 3){
        //     new_node.centerY = node.centerY + node.sideLength * 0.5f; // Lower half
        //     new_node.centerX = node.centerX + node.sideLength * 0.5f; // Left half
        // }
        new_node.centerY = node.centerY  + node.sideLength*(  1.0f - 2.0f*(insert_index<=1));
        new_node.centerX = node.centerX  + node.sideLength*(- 1.0f + 2.0f*(insert_index%2));
        new_node.sideLength = node.sideLength *0.5f; // Halve the side length
        insert_particle(particles, old_particle_idx, tree_array, new_node_idx); // Insert the old particle into the new node
        insert_particle(particles, pidx, tree_array, new_node_idx); // Insert the new particle into the new node
        return;
    }



    
    // int index = current_available_index++;
}
int construct_trees(Node *tree_array, Particle* particles, int pcount, float* bounds){
    current_available_index = PARTICLE_COUNT;
    int tcount = 0;
    for(int i=0; i<pcount; i++){
        float lbx = bounds[0];
        float lby = bounds[1];
        float ubx = bounds[2];
        float uby = bounds[3];
        Particle& particle = particles[i];
        float midX = (lbx + ubx) / 2.0f;
        float midY = (lby + uby) / 2.0f;
        int cell_number = 0;
        // if(particle.x> midX) cell_number += 1;
        // if(particle.y> midY) cell_number += 2;
        cell_number = (particle.x > midX) + 2 * (particle.y > midY);


    }
    

}



void step_particles(Particle* particles, int count){
    // Update positions and velocities
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        particle.x += particle.vx * SUBSTEPS_DT*0.5f;
        particle.y += particle.vy * SUBSTEPS_DT*0.5f;
        particle.vx += particle.ax * SUBSTEPS_DT;
        particle.vy += particle.ay * SUBSTEPS_DT;
        particle.x += particle.vx * SUBSTEPS_DT*0.5f;
        particle.y += particle.vy * SUBSTEPS_DT*0.5f;
        particle.ax = 0.0f; 
        particle.ay = 0.0f;
    }
}

void step_gravity(Particle* particles, int count){
    // Verlet integration for gravity between particles
    for(int i=0; i<count; i++){
        Particle& particle = particles[i];
        float fx = 0.0f;
        float fy = 0.0f;

        for(int j=i; j<count; j++){
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
            other.ax -= force * (dx / dist)/other.mass;
            other.ay -= force * (dy / dist)/other.mass;

        }

    }
    // Update positions and velocities
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
    init_particles(particles, pshape, PARTICLE_COUNT);


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
            // step_gravity(particles, PARTICLE_COUNT);
            compute_forces_naive(particles, PARTICLE_COUNT);
            step_particles(particles, PARTICLE_COUNT);
        }
        sync_shapes(particles, pshape, PARTICLE_COUNT);


        for(int i=0; i<PARTICLE_COUNT; i++){
            window.draw(pshape[i]);
        }

        window.display();
    }
    return 0;
}