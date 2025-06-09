#include <iostream>
#include <stdlib.h>
#include <random>
#include <SFML/Graphics.hpp>

// Simulation properties
#define PARTICLE_COUNT 10000
#define RAD 1.0f
#define dRad 0.0f
#define XMin 000.0f
#define XMax 800.0f
#define YMin 000.0f
#define YMax 800.0f
#define MAGNIFICATION 1.0f
#define MASS 1.0f
#define CENTER_MASS 10000.0f
#define DAMPING 0.4f
#define SUBSTEPS 10



// Numerical properties
#define SUBSTEPS_DT 0.01f
#define MAXIMUM_TREE_SIZE (PARTICLE_COUNT)
#define THETA 1.5f
#define SOFTENING_EPSILON 0.01f 
#define NORM2(x, y) ((x) * (x) + (y) * (y)) 


// Physics properties
#define G 0.10f
#define RAD_MASS_CONSTANT 10.0f
#define MAX_FORCE 10.0f
#define MAX_ACCELERATION 10.0f


// #define DEBUG
// #define DEBUG2
// // #define DEBUG_INSERTION
// #define DEBUG_COM

int current_available_index = PARTICLE_COUNT;
int tree_size = 0;

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
    float com_x;
    float com_y;
    float mass;

#ifdef DEBUG
    Node() : valid(false),
             com_calculated(false),
             leaves{-1, -1, -1, -1}, centerX(0.0f), centerY(0.0f), sideLength(0.0f), com_x(-42.0f), com_y(-42.0f), mass(0.0f) {}
#else
    Node() : valid(false), leaves{-1, -1, -1, -1}, centerX(0.0f), centerY(0.0f), sideLength(0.0f), com_x(-42.0f), com_y(-42.0f), mass(0.0f) {}
#endif
};

Particle *particles_main = new Particle[PARTICLE_COUNT];
sf::CircleShape *pshape_main = new sf::CircleShape[PARTICLE_COUNT];
Node *tree_array_main = new Node[MAXIMUM_TREE_SIZE];
// Helper functions
float randf(float lb = 0.0f, float ub = 1.0f)
{
    return lb + (ub - lb) * (((float)rand()) / ((float)RAND_MAX));
}
void print(auto message)
{
#ifdef DEBUG
    std::cout << message << std::endl;
    #endif
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
        const float radius = MAGNIFICATION * randf( 0.2f * (XMax - XMin) *0.5f, 0.7f * (XMax - XMin) *0.5f); // Use half the width as radius
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

void get_bounding_box(const Particle *particles, int count, float *bounds)
{
    print("Getting bounding box for particles...");
    float minX = particles[0].x;
    float minY = particles[0].y;
    float maxX = particles[0].x;
    float maxY = particles[0].y;

    for (int i = 1; i < count; i++)
    {
        const Particle &particle = particles[i];
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
    bounds[2] = ((maxX - minX) > (maxY - minY) ? (maxX - minX) : (maxY - minY)) ; // Square grid

    // bounds[0] = minX;
    // bounds[1] = minY;
    // bounds[2] = maxX;
    // bounds[3] = maxY;
    print("Bounding box calculaed successfully");
}

void insert_particle(const Particle *particles, Node *tree_array, const int pidx, const int nidx)
{
    // print("Inserting particle " + std::to_string(pidx) + " into node " + std::to_string(nidx));
    Node &node = tree_array[nidx];
#ifdef DEBU
    if (!node.valid)
        print("Node is not valid, cannot insert particle.");
#endif
    // Particle particle = particles[pidx];
    int insert_index = (particles[pidx].x > node.centerX) + 2 * (particles[pidx].y > node.centerY);
    #ifdef DEBUG_INSERTION
    print("particle " + std::to_string(pidx) + " -> node " + std::to_string(nidx) + " || index " + std::to_string(insert_index) + ", center: (" + std::to_string(node.centerX) + ", " + std::to_string(node.centerY) + "), side length: " + std::to_string(node.sideLength));
    #endif
    if (node.leaves[insert_index] == -1)
    {
        #ifdef DEBUG_INSERTION
        print("Node empty, inserting particle " + std::to_string(pidx) + " -> node " + std::to_string(nidx) + " at index " + std::to_string(insert_index));
        // exit(1);
        print("Particle " + std::to_string(pidx) + " inserted successfully at index " + std::to_string(insert_index));
        print("");
        #endif
        node.leaves[insert_index] = pidx;
        return;
    }
    else if (node.leaves[insert_index] >= PARTICLE_COUNT) // It has a node. Go to the node pointed to
    {
        #ifdef DEBUG_INSERTION
        print("Node is not empty, inserting particle " + std::to_string(pidx) + " into child node " + std::to_string(node.leaves[insert_index] - PARTICLE_COUNT));
        #endif
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
            print("Current available index: " + std::to_string(current_available_index));
            exit(1);
        }
#endif
        node.leaves[insert_index] = new_node_idx + PARTICLE_COUNT; // Update the current node to point to the new node
        Node &new_node = tree_array[new_node_idx];
        new_node.valid = true;
        new_node.centerY = node.centerY + node.sideLength * 0.25f * (1 - 2 * (insert_index <= 1));
        new_node.centerX = node.centerX + node.sideLength * 0.25f * (2 * (insert_index % 2) - 1);
        new_node.sideLength = node.sideLength * 0.5f;                           // Halve the side length
        insert_particle(particles, tree_array, old_particle_idx, new_node_idx); // Insert the old particle into the new node
        insert_particle(particles, tree_array, pidx, new_node_idx);             // Insert the new particle into the new node
        return;
    }

    // int index = current_available_index++;
}
Node* construct_trees(const Particle *particles, Node *tree_array, const int pcount, float *bounds)
{
    print("Deleting old tree array...");
    delete[] tree_array;
    print("Constructing tree inner...");
    tree_array = new Node[MAXIMUM_TREE_SIZE];
    print("Tree array constructed successfully. Size: " + std::to_string(MAXIMUM_TREE_SIZE));
    tree_array[0].valid = true;
    tree_array[0].centerX = bounds[0];
    tree_array[0].centerY = bounds[1];
    tree_array[0].sideLength = bounds[2];
    tree_array[0].leaves[0] = -1;
    tree_array[0].leaves[1] = -1;
    tree_array[0].leaves[2] = -1;
    tree_array[0].leaves[3] = -1;
    
    current_available_index = 1;
    print("Inserting particles into tree...");
    for (int i = 0; i < pcount; i++)
    {
        insert_particle(particles, tree_array, i, 0);
    }
    tree_size = current_available_index;
    print("Particles inserted successfully. Current available index: " + std::to_string(tree_size));
    print("testing tree_array[0].valid: " + std::to_string(tree_array[tree_size-1].valid));
    return tree_array;
}

void computer_tree_coms(Particle *particles, Node *tree_array, const int pcount, const int ncount)
{
    print(tree_array[0].valid ? "Root node is valid." : "Root node is not valid.");
    print("Entering inner com function \n");
    for (int i = ncount - 1; i >= 0; i--)
    {

        #ifdef DEBUG_COM
        if(i>=tree_size)
        {
            print("Node index out of bounds: " + std::to_string(i) + ", tree size: " + std::to_string(tree_size));
            exit(1);
        }

        #endif
        // print("Attempting com compute");
        // print("Computing COM for node");
        Node &node = tree_array[i];

        // print(node.valid ? "Node is valid." : "Node is not valid.");
        // print("Able to access node");
#ifdef DEBUG
        if (!node.valid)
        {
            print("Node is not valid, cannot compute COM.");
            exit(1);
            // continue;
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
                #ifdef DEBUG_COM
                if (leaf - PARTICLE_COUNT >= ncount) print("Leaf index out of bounds: " + std::to_string(leaf - PARTICLE_COUNT));
                // print("Trying to compute COM for child node " + std::to_string(leaf - PARTICLE_COUNT) + " in parent node " + std::to_string(i));
                #endif
                const Node &child_node = tree_array[leaf - PARTICLE_COUNT];
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
                comx += child_node.com_x * child_node.mass;
                comy += child_node.com_y * child_node.mass;
                mass += child_node.mass;
            }
            else if (leaf >= 0)
            { // It has a particle
                const Particle &particle = particles[leaf];
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
#ifdef DEBUG
        if (mass == 0.0f)
        {
            print("Mass is zero, cannot compute COM.");
            exit(1);
        }
#endif
        node.com_x = comx / mass;
        node.com_y = comy / mass;
        node.mass = mass;
#ifdef DEBUG
        node.com_calculated = true;
#endif
    }
}

// void compute_forces_at_point_by_node(const Node* tree_array, const int ncount, )
void compute_accelerations_at_point_bh(const Node *tree_array, const Particle *particles, const int ncount, const int nidx, const float &x, const float &y, const float &theta, float &ax, float &ay)
{
    int index = 0;
    bool leaf = false;
#ifdef DEBUG
    if (!tree_array[nidx].valid)
    {
        print("Node is not valid, cannot compute forces.");
        // exit(1);
        // continue;
    }
#endif
    const float dx = tree_array[nidx].com_x - x;
    const float dy = tree_array[nidx].com_y - y;
    const float inv_dist_sq = 1.0f/(NORM2(dx, dy)  + SOFTENING_EPSILON);
    const float self_theta_sq = (tree_array[nidx].sideLength * tree_array[nidx].sideLength) *inv_dist_sq;
    if (self_theta_sq < theta * theta)
    {
        // If the node is far enough, treat it as a single particle
        const float inv_dist = sqrtf(inv_dist_sq);
#ifdef DEBUG2
        if (inv_dist > 1e2f)
            return;
#endif
        const float acc = G * (tree_array[nidx].mass) * inv_dist_sq;
        ax += acc * (dx * inv_dist);
        ay += acc * (dy * inv_dist);
        return;
    }
    else
    {
        for (int i = 0; i < 4; i++)
        {
            int child_index = tree_array[nidx].leaves[i];
       
            
            if (child_index == -1)
                continue;
            if (child_index >= PARTICLE_COUNT)
            {
                #ifdef DEBUG
                if(child_index - PARTICLE_COUNT >= ncount)
                {
                    print("Child index out of bounds: " + std::to_string(child_index - PARTICLE_COUNT));
                    exit(1);
                }
                #endif
                compute_accelerations_at_point_bh(tree_array, particles, ncount, child_index - PARTICLE_COUNT, x, y, theta, ax, ay);
            }
            else if (child_index >= 0)
            {
                // It's a particle
                const Particle &particle = particles[child_index];
                const float dx_p = particle.x - x;
                const float dy_p = particle.y - y;
                const float dist_sq_p = NORM2(dx_p, dy_p)+ SOFTENING_EPSILON;
                if (dist_sq_p < 1e-2f)
                    continue; // Avoid division by zero
                const float inv_dist_p = 1.0f/sqrtf(dist_sq_p);
                const float acc_p = G * (particle.mass) / dist_sq_p;
                ax += acc_p * (dx_p * inv_dist_p);
                ay += acc_p * (dy_p * inv_dist_p);
            } 
            else
            {
                print("Unknown case in compute_accelerations_at_point_bh, child index: " + std::to_string(child_index));
                exit(1);
            }
        }
    }
}
void calculate_accelerations_barnes_hut(Particle *particles, const Node *tree_array, const int pcount, const int ncount, const float theta)
{
    for (int i = 0; i < pcount; i++)
    {
        Particle &target = particles[i];
        const float x = target.x;
        const float y = target.y;
        float ax = 0.0f;
        float ay = 0.0f;
        compute_accelerations_at_point_bh(tree_array, particles, ncount, 0, x, y, theta, ax, ay);
        target.ax = ax;
        target.ay = ay;
    }
}

void perform_barnes_hut(Particle *particles, Node *tree_array, const int pcount, const float theta){
    float bounds[3];
    print("Performing Barnes-Hut algorithm...");
    get_bounding_box(particles, pcount, bounds);
    print("Constructing trees...");
    tree_array_main = construct_trees(particles, tree_array_main, pcount, bounds);
    print("Tree checking...");
    print(tree_array_main[0].valid ? "Root node is valid." : "Root node is not valid.");
    const int ncount = current_available_index;
    print("Trees constructed successfully. Number of nodes: " + std::to_string(ncount));
    print("Computing center of mass for nodes...");
    computer_tree_coms(particles, tree_array_main, pcount, ncount);
    print("Center of mass computed successfully.");
    print("Calculating accelerations using Barnes-Hut algorithm...");
    calculate_accelerations_barnes_hut(particles, tree_array_main, pcount, ncount, theta);
    print("Accelerations calculated successfully.");
}

void step_particles(Particle *particles, int count)
{
    // Update positions and velocities
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        float mag_acc_sq = NORM2(particle.ax, particle.ay);
        if (mag_acc_sq > MAX_ACCELERATION * MAX_ACCELERATION){
            // Limit acceleration to MAX_ACCELERATION
            float mag_acc_inv = 1/sqrtf(mag_acc_sq);
            particle.ax = (particle.ax * mag_acc_inv) * MAX_ACCELERATION;
            particle.ay = (particle.ay * mag_acc_inv) * MAX_ACCELERATION;
        }


        particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
        particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
        particle.vx += particle.ax * SUBSTEPS_DT;
        particle.vy += particle.ay * SUBSTEPS_DT;
        particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
        particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
        particle.ax = 0.0f;
        particle.ay = 0.0f;



        if(particle.x < XMin*MAGNIFICATION) {
            particle.x = XMin * MAGNIFICATION;
            particle.vx = -DAMPING*particle.vx;
        }
        if(particle.x > XMax*MAGNIFICATION) {
            particle.x = XMax * MAGNIFICATION;
            particle.vx = -DAMPING*particle.vx;
        }
        if(particle.y < YMin*MAGNIFICATION) {
            particle.y = YMin * MAGNIFICATION;
            particle.vy = -DAMPING*particle.vy;
        }
        if(particle.y > YMax*MAGNIFICATION) {
            particle.y = YMax * MAGNIFICATION;
            particle.vy = -DAMPING*particle.vy;
        }
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
    
    init_particles(particles_main, pshape_main, PARTICLE_COUNT);
    print("Particles initialized successfully. Particle count: " + std::to_string(PARTICLE_COUNT));
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
            // compute_forces_naive(particles_main, PARTICLE_COUNT);
            print("Performing Barnes-Hut algorithm...");
            perform_barnes_hut(particles_main, tree_array_main, PARTICLE_COUNT, THETA);
            step_particles(particles_main, PARTICLE_COUNT);

        }
        sync_shapes(particles_main, pshape_main, PARTICLE_COUNT);

        for (int i = 0; i < PARTICLE_COUNT; i++)
        {
            window.draw(pshape_main[i]);
        }

        window.display();
    }
    return 0;
}