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

#define MAXIMUM_TREE_SIZE (5 * PARTICLE_COUNT)
// Physics properties
#define G 0.10f
#define RAD_MASS_CONSTANT 10.0f
#define MAX_FORCE 1000.0f

#define DEBUG
#define DEBUG2

int current_available_index = PARTICLE_COUNT;

// //Memory objects

struct Particle
{
    float x;
    float y;
    float vx;
    float vy;
    float ax;
    float ay;
    float rad;
    float mass;
};

struct Node
{

    bool valid;
#ifdef DEBUG
    bool com_calculated;
#endif
    int leaves[4];
    float centerX;
    float centerY;
    float sideLength;
    float mass_x;
    float mass_y;
    float mass;

#ifdef DEBUG
    Node() : valid(false),
             com_calculated(false),
             leaves{-1, -1, -1, -1}, centerX(0.0f), centerY(0.0f), sideLength(0.0f), mass_x(-42.0f), mass_y(-42.0f), mass(0.0f) {}
#else
    Node() : valid(false), leaves{-1, -1, -1, -1}, centerX(0.0f), centerY(0.0f), sideLength(0.0f), mass_x(-42.0f), mass_y(-42.0f), mass(0.0f) {}
#endif
};

Particle *particles = new Particle[PARTICLE_COUNT];
sf::CircleShape *pshape = new sf::CircleShape[PARTICLE_COUNT];

// Helper functions
float randf(float lb = 0.0f, float ub = 1.0f)
{
    return lb + (ub - lb) * (((float)rand()) / ((float)RAND_MAX));
}
void print(auto message)
{
    std::cout << message << std::endl;
}

// Supplementary functions

void init_particles(Particle *particles, sf::CircleShape *shapes, int count)
{
    const float centerX = (XMax + XMin) / 2.0f;
    const float centerY = (YMax + YMin) / 2.0f;
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        sf::CircleShape &shape = shapes[i];
        // particle.x = randf(XMin*MAGNIFICATION, XMax*MAGNIFICATION);
        // particle.y = randf(YMin*MAGNIFICATION, YMax*MAGNIFICATION);
        // Initialize velocity to be in stable circular motion. centripital == gravitational force
        float angle = randf(0.0f, 2.0f * 3.14159265358979323846f);
        const float radius = MAGNIFICATION * randf(0.0f, 0.5f * (XMax - XMin) / 2.0f); // Use half the width as radius
        particle.x = MAGNIFICATION * ((XMax + XMin) / 2.0f + radius * cosf(angle));
        particle.y = MAGNIFICATION * ((YMax + YMin) / 2.0f + radius * sinf(angle));
        // particle.x = MAGNIFICATION * (XMax + XMin) / 2.0f + RAD * cosf(angle);
        // particle.y = MAGNIFICATION * (YMax + YMin) / 2.0f + RAD * sinf(angle);
        float speed = sqrtf(G * CENTER_MASS / (radius)); // Speed for stable circular motion
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
    particles[0].x = MAGNIFICATION * (XMax + XMin) / 2.0f;
    particles[0].y = MAGNIFICATION * (YMax + YMin) / 2.0f;
    particles[0].rad = 3 * RAD;
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
void compute_forces_naive(Particle *particles, int count)
{
    // Compute forces between particles
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        float fx = 0.0f;
        float fy = 0.0f;

        for (int j = 0; j < count; j++)
        {
            if (i == j)
                continue; // Skip self

            Particle &other = particles[j];
            float dx = (other.x - particle.x);
            float dy = (other.y - particle.y);
            float dist_sq = dx * dx + dy * dy;

            if (dist_sq < 1e-2f)
                continue; // Avoid division by zero

            const float dist = sqrtf(dist_sq);
            float force = G * (particle.mass * other.mass) / dist_sq;
            force = force > MAX_FORCE ? 0 : force; // Limit force to MAX_FORCE
            particle.ax += force * (dx / dist) / particle.mass;
            particle.ay += force * (dy / dist) / particle.mass;
        }
    }
}

// void insert_point(int *tree_array, int)

void get_bounding_box(Particle *particles, int count, float *bounds)
{
    float minX = particles[0].x;
    float minY = particles[0].y;
    float maxX = particles[0].x;
    float maxY = particles[0].y;

    for (int i = 1; i < count; i++)
    {
        Particle &particle = particles[i];
        if (particle.x < minX)
            minX = particle.x;
        if (particle.y < minY)
            minY = particle.y;
        if (particle.x > maxX)
            maxX = particle.x;
        if (particle.y > maxY)
            maxY = particle.y;
    }
    bounds[0] = (minX + maxX) * 0.5f;
    bounds[1] = (minY + maxY) * 0.5f;
    bounds[2] = ((maxX - minX) > (maxY - minY) ? (maxX - minX) : (maxY - minY)) * 0.5f; // Square grid

    // bounds[0] = minX;
    // bounds[1] = minY;
    // bounds[2] = maxX;
    // bounds[3] = maxY;
}

void insert_particle(const Particle *particles, Node *tree_array, const int pidx, const int index)
{
    Node &node = tree_array[index];
#ifdef DEBUG
    if (!node.valid)
        print("Node is not valid, cannot insert particle.");
#endif
    // Particle particle = particles[pidx];
    int insert_index = (particles[pidx].x > node.centerX) + 2 * (particles[pidx].y > node.centerY);
    if (node.leaves[insert_index] == -1)
    {
        node.leaves[insert_index] = pidx;
        return;
    }
    else if (node.leaves[insert_index] >= PARTICLE_COUNT) // It has a node. Go to the node pointed to
    {
        insert_particle(particles, tree_array, pidx, node.leaves[insert_index] - PARTICLE_COUNT);
        return;
    }

    else if (node.leaves[insert_index] < PARTICLE_COUNT)
    {
        // Node is a leaf, we need to split it
        int old_particle_idx = node.leaves[insert_index];
        int new_node_idx = current_available_index++;
#ifdef DEBUG
        if (current_available_index >= MAXIMUM_TREE_SIZE)
        {
            print("Maximum tree size reached, cannot insert more particles.");
            exit(1);
        }
#endif
        node.leaves[insert_index] = new_node_idx + PARTICLE_COUNT; // Update the current node to point to the new node
        Node &new_node = tree_array[new_node_idx];
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
        new_node.centerY = node.centerY + node.sideLength * (1.0f - 2.0f * (insert_index <= 1));
        new_node.centerX = node.centerX + node.sideLength * (-1.0f + 2.0f * (insert_index % 2));
        new_node.sideLength = node.sideLength * 0.5f;                           // Halve the side length
        insert_particle(particles, tree_array, old_particle_idx, new_node_idx); // Insert the old particle into the new node
        insert_particle(particles, tree_array, pidx, new_node_idx);             // Insert the new particle into the new node
        return;
    }

    // int index = current_available_index++;
}
int construct_trees(Node *tree_array, const Particle *particles, const int pcount, float *bounds)
{
    delete[] tree_array;
    tree_array = new Node[MAXIMUM_TREE_SIZE];
    current_available_index = 0;
    tree_array[0].valid = true;
    tree_array[0].centerX = bounds[0];
    tree_array[0].centerY = bounds[1];
    tree_array[0].sideLength = bounds[2];
    tree_array[0].leaves[0] = -1;
    tree_array[0].leaves[1] = -1;
    tree_array[0].leaves[2] = -1;
    tree_array[0].leaves[3] = -1;

    for (int i = 0; i < pcount; i++)
    {
        insert_particle(particles, tree_array, i, 0);
    }
    return current_available_index;
}

void computer_tree_coms(const Particle *particle, Node *tree_array, const int pcount, const int ncount)
{
    for (int i = ncount - 1; i >= 0; i--)
    {
        Node &node = tree_array[i];
#ifdef DEBUG
        if (!node.valid)
        {
            print("Node is not valid, cannot compute COM.");
            exit(1);
        }
#endif
        float comx = 0.0f;
        float comy = 0.0f;
        float mass = 0.0f;
        for (int j = 0; j < 4; j++)
        {
            const int &leaf = node.leaves[j];
            if (leaf == -1)
                continue;
            if (leaf >= PARTICLE_COUNT)
            { // It has a node
                Node &child_node = tree_array[leaf - PARTICLE_COUNT];
#ifdef DEBUG
                if (!child_node.valid)
                {
                    print("Child node is not valid, cannot compute COM.");
                    exit(1);
                }
                if (!child_node.com_calculated)
                {
                    print("Child node COM not calculated, cannot compute parent COM.");
                    exit(1);
                }
#endif
                comx += child_node.mass_x * child_node.mass;
                comy += child_node.mass_y * child_node.mass;
                mass += child_node.mass;
            }
            else if (leaf >= 0)
            { // It has a particle
                Particle &particle = particles[node.leaves[j]];
                comx += particle.x * particle.mass;
                comy += particle.y * particle.mass;
                mass += particle.mass;
            }
            else
            {
                print("Unknown case in computer_tree_coms, leaf index: " + std::to_string(leaf));
                exit(1);
            }
        }
        node.mass_x = comx / mass;
        node.mass_y = comy / mass;
        node.mass = mass;
#ifdef DEBUG
        node.com_calculated = true;
#endif
    }
}

// void compute_forces_at_point_by_node(const Node* tree_array, const int ncount, )
void compute_forces_at_point_bh(const Node *tree_array, const int ncount, const int nidx, const float &x, const float &y, const float &theta, float &fx, float &fy)
{
    int index = 0;
    bool leaf = false;
#ifdef DEBUG
    if (!tree_array[nidx].valid)
    {
        print("Node is not valid, cannot compute forces.");
        exit(1);
    }
#endif
    const float dx = tree_array[nidx].mass_x - x;
    const float dy = tree_array[nidx].mass_y - y;
    const float dist_sq = dx * dx + dy * dy;
    const float self_theta_sq = (tree_array[nidx].sideLength*tree_array[nidx].sideLength)/dist_sq;
    if (self_theta_sq < theta * theta){
        // If the node is far enough, treat it as a single particle
        const float dist = sqrtf(dist_sq);
        #ifdef DEBUG2
        if (dist < 1e-2f) return;
        #endif
        const float force = G * (tree_array[nidx].mass) / dist_sq;
        fx += force * (dx / dist);
        fy += force * (dy / dist);
        return;
    }else{
        for(int i=0; i<4; i++){
            int child_index = tree_array[nidx].leaves[i];
            if(child_index == -1) continue; 
            if(child_index >= PARTICLE_COUNT){
                compute_forces_at_point_bh(tree_array, ncount, child_index - PARTICLE_COUNT, x, y, theta, fx, fy);
            }else if(child_index >=0){
                // It's a particle
                const Particle &particle = particles[child_index];
                const float dx_p = particle.x - x;
                const float dy_p = particle.y - y;
                const float dist_sq_p = dx_p * dx_p + dy_p * dy_p;
                if (dist_sq_p < 1e-2f) continue; // Avoid division by zero
                const float dist_p = sqrtf(dist_sq_p);
                const float force_p = G * (particle.mass) / dist_sq_p;
                fx += force_p * (dx_p / dist_p);
                fy += force_p * (dy_p / dist_p);
            }
            else
            {
                print("Unknown case in compute_forces_at_point_bh, child index: " + std::to_string(child_index));
                exit(1);
            }
        }
    }

}
void calculate_forces_barnes_hut(Particle *particles, const Node *tree_array, const int pcount, const int ncount, const float theta)
{
    for (int i = 0; i < pcount; i++)
    {
        Particle &target = particles[i];
        const float x = target.x;
        const float y = target.y;
        // Todo
    }
}

void step_particles(Particle *particles, int count)
{
    // Update positions and velocities
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
        particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
        particle.vx += particle.ax * SUBSTEPS_DT;
        particle.vy += particle.ay * SUBSTEPS_DT;
        particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
        particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
        particle.ax = 0.0f;
        particle.ay = 0.0f;
    }
}

void step_gravity(Particle *particles, int count)
{
    // Verlet integration for gravity between particles
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        float fx = 0.0f;
        float fy = 0.0f;

        for (int j = i; j < count; j++)
        {
            if (i == j)
                continue; // Skip self

            Particle &other = particles[j];
            float dx = (other.x - particle.x);
            float dy = (other.y - particle.y);
            float dist_sq = dx * dx + dy * dy;

            if (dist_sq < 1e-2f)
                continue; // Avoid division by zero

            const float dist = sqrtf(dist_sq);
            float force = G * (particle.mass * other.mass) / dist_sq;
            force = force > MAX_FORCE ? 0 : force; // Limit force to MAX_FORCE
            particle.ax += force * (dx / dist) / particle.mass;
            particle.ay += force * (dy / dist) / particle.mass;
            other.ax -= force * (dx / dist) / other.mass;
            other.ay -= force * (dy / dist) / other.mass;
        }
    }
    // Update positions and velocities
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
        particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
        // if(i==100) std::cout << "Particle " << i << " position: (" << particle.x/MAGNIFICATION << ", " << particle.y/MAGNIFICATION << ")" << std::endl;
        particle.vx += particle.ax * SUBSTEPS_DT;
        particle.vy += particle.ay * SUBSTEPS_DT;
        particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
        particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
        particle.ax = 0.0f;
        particle.ay = 0.0f;
    }
}

void sync_shapes(Particle *particles, sf::CircleShape *shapes, int count)
{
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        sf::CircleShape &shape = shapes[i];
        shape.setPosition(particle.x / MAGNIFICATION, particle.y / MAGNIFICATION);
    }
}

int main()
{
    // Initialize random seed
    srand(42);
    init_particles(particles, pshape, PARTICLE_COUNT);

    sf::RenderWindow window(sf::VideoMode(XMax - XMin, YMax - YMin), "Particle Simulation");
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
        for (int i = 0; i < SUBSTEPS; i++)
        {
            // Perform substeps
            // step_gravity(particles, PARTICLE_COUNT);
            compute_forces_naive(particles, PARTICLE_COUNT);
            step_particles(particles, PARTICLE_COUNT);
        }
        sync_shapes(particles, pshape, PARTICLE_COUNT);

        for (int i = 0; i < PARTICLE_COUNT; i++)
        {
            window.draw(pshape[i]);
        }

        window.display();
    }
    return 0;
}