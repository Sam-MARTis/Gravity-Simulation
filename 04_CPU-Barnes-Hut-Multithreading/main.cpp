#include <iostream>
#include <stdlib.h>
#include <random>
#include <SFML/Graphics.hpp>
#include <omp.h>
#include <chrono>

// g++ main.cpp -o grav -lsfml-graphics -lsfml-window -lsfml-system -fopenmp

// Simulation properties
#define PARTICLE_COUNT 10000
#define CLOUD_VEL_MULTIPLIER 1.0f
#define RAD 0.1f
#define dRad 0.05f
#define XMin 000.0f
#define XMax 900.0f
#define YMin 000.0f
#define YMax 900.0f
#define INNER_RADII 0.05f
#define OUTER_RADII 0.7f
#define COLLISION_RADIUS_MULTIPLER 2.0f

#define MAGNIFICATION 1.0f
#define MASS 100.0f
#define CENTER_MASS 3000000.0f
#define CENTER_RAD 2.0f
#define CENTER_MASS_COUNT 2
#define D_CENTER 3.0f
#define DAMPING 0.4f
#define SUBSTEPS 5

#define INTERNAL_ELASTICITY 1.0

// Numerical properties
#define SUBSTEPS_DT 0.001
#define MAXIMUM_TREE_SIZE (5 * PARTICLE_COUNT)
#define THETA 1.5f
#define SOFTENING_EPSILON 0.01f
#define NORM2(x, y) ((x) * (x) + (y) * (y))
#define SQ(x) ((x) * (x))
#define MAX_FORCE 10.0f
#define MAX_VEL_IMPULSE 40.0f
#define MAX_ACCELERATION 5000000.0f

#define COLLISION_HANDLING_PER_ITERATION 10
#define MAX_COLLISION_SEPRATION_FRACTION 0.3

#define TIMEIT(func, message) \
    {                         \
        func;                 \
    }

// Timer timer = Timer(message, 0);       \ 

// Physics properties
#define G (1.0f)
#define RAD_MASS_CONSTANT 10.0f
#define RAD_LUMINOSITY 4050.0f

// #define DEBUG
// #define DEBUG2
// // #define DEBUG_INSERTION
// #define DEBUG_COM

int current_available_index = PARTICLE_COUNT;
int tree_size = 0;

// //Memory objects

struct Particle
{
    double x;
    double y;
    double x_prev;
    double y_prev;
    // double vx;
    // double vy;
    double ax;
    double ay;
    float rad;
    float mass;
    double dx;
    double dy;
    // double dvx;
    // double dvy;
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

struct Timer
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<float> duration;
    std::string name;
    int num;

    Timer(std::string in, int ran = 1)
    {
        start = std::chrono::high_resolution_clock::now();
        name = in;
        num = ran;
    }
    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        if (num % 10 == 0)
        {
            std::cout << name << " took: " << duration.count() * 1000.0f << "ms\n";
        }
    }
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

void print2(auto message)
{
    std::cout << message << std::endl;
}

// Supplementary functions

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
    bounds[2] = ((maxX - minX) > (maxY - minY) ? (maxX - minX) : (maxY - minY)); // Square grid

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
Node *construct_trees(const Particle *particles, Node *tree_array, const int pcount, float *bounds)
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
    print("testing tree_array[0].valid: " + std::to_string(tree_array[tree_size - 1].valid));
    return tree_array;
}

void computer_tree_coms(Particle *particles, Node *tree_array, const int pcount, const int ncount)
{
    print(tree_array[0].valid ? "Root node is valid." : "Root node is not valid.");
    print("Entering inner com function \n");
    // #pragma omp parallel for
    for (int i = ncount - 1; i >= 0; i--)
    {

#ifdef DEBUG_COM
        if (i >= tree_size)
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
                if (leaf - PARTICLE_COUNT >= ncount)
                    print("Leaf index out of bounds: " + std::to_string(leaf - PARTICLE_COUNT));
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
void compute_accelerations_at_particle_bh(const Node *tree_array, Particle *particles, const int ncount, const int nidx, const int &pidx, const float &theta, double &ax, double &ay)
{
    int index = 0;
    bool leaf = false;
    Particle &particle_original = particles[pidx];
    const double x = particle_original.x;
    const double y = particle_original.y;

#ifdef DEBUG
    if (!tree_array[nidx].valid)
    {
        print("Node is not valid, cannot compute forces.");
        // exit(1);
        // continue;
    }
#endif
    const double dx = tree_array[nidx].com_x - x;
    const double dy = tree_array[nidx].com_y - y;
    const double inv_dist_sq = 1.0f / (NORM2(dx, dy) + SOFTENING_EPSILON);
    const float self_theta_sq = (tree_array[nidx].sideLength * tree_array[nidx].sideLength) * inv_dist_sq;
    if (self_theta_sq < theta * theta)
    {
        // If the node is far enough, treat it as a single particle
        const float inv_dist = sqrtf(inv_dist_sq);
#ifdef DEBUG2
        if (inv_dist > 1e2f)
            return;
#endif
        const double acc = G * (tree_array[nidx].mass) * inv_dist_sq;
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
                if (child_index - PARTICLE_COUNT >= ncount)
                {
                    print("Child index out of bounds: " + std::to_string(child_index - PARTICLE_COUNT));
                    exit(1);
                }
#endif
                compute_accelerations_at_particle_bh(tree_array, particles, ncount, child_index - PARTICLE_COUNT, pidx, theta, ax, ay);
            }
            else if (child_index >= 0)
            {
                // It's a particle
                Particle &particle_other = particles[child_index];
                const double dx_p = particle_other.x - x;
                const double dy_p = particle_other.y - y;
                const double dist_sq_p = NORM2(dx_p, dy_p) + SOFTENING_EPSILON;

                // if (dist_sq_p < 1e-2f)
                //     continue; // Avoid division by zero
                const double inv_dist_p = 1.0f / sqrtf(dist_sq_p);
                // if ((dist_sq_p - SOFTENING_EPSILON) < SQ(particle_other.rad + particle_original.rad))
                // {
                //     // // Handle collision
                //     // const float relavtive_distance_inv = sqrtf(dist_sq_p);
                //     const double dR = (particle_other.rad + particle_original.rad) - 1 / inv_dist_p;
                //     const double mass_ratio = particle_other.mass / (particle_other.mass + particle_original.mass);
                //     const double relative_vx = particle_other.vx - particle_original.vx;
                //     const double relative_vy = particle_other.vy - particle_original.vy;
                //     const double dV = ((relative_vx * dx_p + relative_vy * dy_p) * inv_dist_p) * (1 + INTERNAL_ELASTICITY) * (mass_ratio);

                //     particle_original.dx -= (dx_p * inv_dist_p) * dR * mass_ratio;
                //     particle_original.dy -= (dy_p * inv_dist_p) * dR * mass_ratio;

                //     particle_original.dvx += (dx_p * inv_dist_p) * dV;
                //     particle_original.dvy += (dy_p * inv_dist_p) * dV;
                // }
                // else
                // {
                const double acc_p = G * (particle_other.mass) * SQ(inv_dist_p);
                ax += acc_p * (dx_p * inv_dist_p);
                ay += acc_p * (dy_p * inv_dist_p);
                // }
            }
            else
            {
#ifdef DEBUG
                print("Unknown case in compute_accelerations_at_point_bh, child index: " + std::to_string(child_index));
#endif
                exit(1);
            }
        }
    }
}
void calculate_accelerations_barnes_hut(Particle *particles, const Node *tree_array, const int pcount, const int ncount, const float theta)
{
#pragma omp parallel for
    for (int i = 0; i < pcount; i++)
    {
        Particle &target = particles[i];
        const double x = target.x;
        const double y = target.y;
        double ax = 0.0;
        double ay = 0.0;
        compute_accelerations_at_particle_bh(tree_array, particles, ncount, 0, i, theta, ax, ay);
        target.ax = ax;
        target.ay = ay;
    }
}
void apply_particle_impulses(Particle *particles, const int pcount)
{
    // Apply impulses to particles
    for (int i = 0; i < pcount; i++)
    {
        Particle &particle = particles[i];
        particle.x += particle.dx;
        particle.y += particle.dy;
        // particle.vx += particle.dvx;
        // particle.vy += particle.dvy;
        particle.dx = 0.0f;
        particle.dy = 0.0f;
        // particle.dvx = 0.0f;
        // particle.dvy = 0.0f;
    }
}
void handle_collision_particle(Particle *particles, const Node *tree_array, const int pcount, const int ncount, const int pidx, const int nidx)
{
    Particle &particle = particles[pidx];
    const Node &node = tree_array[nidx];
    // print2("A");
    const bool left = (particle.x - particle.rad * COLLISION_RADIUS_MULTIPLER) < node.centerX;
    const bool top = (particle.y - particle.rad * COLLISION_RADIUS_MULTIPLER) < node.centerY;
    const bool right = (particle.x + particle.rad * COLLISION_RADIUS_MULTIPLER) > node.centerX;
    const bool bottom = (particle.y + particle.rad * COLLISION_RADIUS_MULTIPLER) > node.centerY;
    // print2("B");
    bool is_intersecting[4] = {left && top, right && top, left && bottom, right && bottom};
#ifdef DEBUG
    if (!node.valid)
    {
        print("Node is not valid, cannot handle collisions.");
        exit(1);
    }
#endif
    for (int i = 0; i < 4; i++)
    {
        const int &leaf = node.leaves[i];
        if (leaf == -1)
            continue;
        if (!is_intersecting[i])
            continue;
        if (leaf >= PARTICLE_COUNT)
        {
            // Node
            if (leaf - PARTICLE_COUNT >= ncount)
            {
                print("Leaf index out of bounds: " + std::to_string(leaf - PARTICLE_COUNT));
                exit(1);
            }
            handle_collision_particle(particles, tree_array, pcount, ncount, pidx, leaf - PARTICLE_COUNT);
            continue;
        }
        if (leaf >= 0)
        {
            // Particle
            // print2("C");

            /*
            void handle_collision(float *pos, const int &id1, const int &id2)
                {
                    const int p1 = (id1) * 2;
                    const int p2 = (id2) * 2;
                    const float dx = pos[p2] - pos[p1];
                    const float dy = pos[p2 + 1] - pos[p1 + 1];
                    const float distance_squared = dx*dx + dy*dy;
                    if (distance_squared > RADIUS_SQUARED_TIMES_FOUR)
                    {
                        return;
                    }
                    const float distance = sqrtf(distance_squared);
                    if (distance < 1e-6f)
                        return;
                    const float nx = dx / distance;
                    const float ny = dy / distance;
                    const float dR = RADIUS - 0.5f * distance;
                    pos[p1] -= (nx * dR);
                    pos[p1 + 1] -= (ny * dR);
                    pos[p2] += (nx * dR);
                    pos[p2 + 1] += (ny * dR);
                }

            */
            Particle &other = particles[leaf];
            if (leaf == pidx)
                continue;
            const double dx = other.x - particle.x;
            const double dy = other.y - particle.y;
            const double dist_sq = NORM2(dx, dy);
            if (dist_sq < SQ(particle.rad + other.rad))
            {
                // Handle collision
                const double inv_dist = 1.0 / sqrtf(dist_sq);
                double dR = (particle.rad + other.rad) - 1.0 / inv_dist;
                const double max_dr_allowed = MAX_COLLISION_SEPRATION_FRACTION * (particle.rad + other.rad);
                dR = dR >max_dr_allowed? max_dr_allowed : dR; // Limit separation to a fraction of the radius

                const double mass_ratio = other.mass / (other.mass + particle.mass);
                // std::cout<<mass_ratio<<std::endl;
                // const double relative_vx = other.vx - particle.vx;
                // const double relative_vy = other.vy - particle.vy;
                // const double dV = ((relative_vx * dx + relative_vy * dy) * inv_dist) * (1 + INTERNAL_ELASTICITY) * (mass_ratio);
                const double nx = dx * inv_dist;
                const double ny = dy * inv_dist;
                particle.x -= (nx) * dR * mass_ratio*(1+INTERNAL_ELASTICITY);
                particle.y -= (ny) * dR * mass_ratio*(1+INTERNAL_ELASTICITY);
                other.x += (nx) * dR * (1 - mass_ratio)*(1+INTERNAL_ELASTICITY);
                other.y += (ny) * dR * (1 - mass_ratio)*(1+INTERNAL_ELASTICITY);
                // particle.dvx += (dx * inv_dist) * dV;
                // particle.dvy += (dy * inv_dist) * dV;
            }
        }
    }
}
void handle_collisions(Particle *particles, const Node *tree_array, const int pcount, const int ncount, const int iterations)
{
#pragma omp parallel 
{

        for (int iter = 0; iter < iterations; iter++){
            #pragma omp for 
                for (int i = 0; i < pcount; i++){
                    Particle &particle = particles[i];
                    const double x = particle.x;
                    const double y = particle.y;
                    const double rad = particle.rad * COLLISION_RADIUS_MULTIPLER;
                    const double rad_sq = rad * rad;
#ifdef DEBUG
                    if (x < tree_array[0].centerX - tree_array[0].sideLength * 0.5f || x > tree_array[0].centerX + tree_array[0].sideLength * 0.5f ||
                        y < tree_array[0].centerY - tree_array[0].sideLength * 0.5f || y > tree_array[0].centerY + tree_array[0].sideLength * 0.5f)
                    {
                        print("Not within root bounds. Shouldnt happen.");
                        exit(1);
                    }
#endif

                    handle_collision_particle(particles, tree_array, pcount, ncount, i, 0);
                    }
                
            #pragma omp for 
                for (int i = 0; i < pcount; i++)
                {
                    Particle &particle = particles[i];
                    particle.x += particle.dx;
                    particle.y += particle.dy;
                    // particle.vx += particle.dvx;
                    // particle.vy += particle.dvy;
                    particle.dx = 0.0f;
                    particle.dy = 0.0f;
                    // particle.dvx = 0.0f;
                    // particle.dvy = 0.0f;
                
                }
            
        
            }
        }
}
void stop_particles(Particle *particles, int count)
{
    // Stop particles by setting their velocities to zero
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        particle.ax = 0.0f;
        particle.ay = 0.0f;
        particle.dx = 0.0f;
        particle.dy = 0.0f;
        particle.x_prev = particle.x;
        particle.y_prev = particle.y;
    }
}

void step_particles(Particle *particles, int count)
{
    // Update positions and velocities
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        float mag_acc_sq = NORM2(particle.ax, particle.ay);
        if (mag_acc_sq > MAX_ACCELERATION * MAX_ACCELERATION)
        {
            // Limit acceleration to MAX_ACCELERATION
            float mag_acc_inv = 1 / sqrtf(mag_acc_sq);
            particle.ax = (particle.ax * mag_acc_inv) * MAX_ACCELERATION;
            particle.ay = (particle.ay * mag_acc_inv) * MAX_ACCELERATION;
        }

        // particle.x += particle.vx * SUBSTEPS_DT * 0.5;
        // particle.y += particle.vy * SUBSTEPS_DT * 0.5;
        // particle.vx += particle.ax * SUBSTEPS_DT;
        // particle.vy += particle.ay * SUBSTEPS_DT;
        // particle.x += particle.vx * SUBSTEPS_DT * 0.5;
        // particle.y += particle.vy * SUBSTEPS_DT * 0.5;
        // Verlet integration
        const float oldx = particle.x;
        const float oldy = particle.y;
        particle.x = 2 * particle.x - particle.x_prev + particle.ax * SUBSTEPS_DT * SUBSTEPS_DT;
        particle.y = 2 * particle.y - particle.y_prev + particle.ay * SUBSTEPS_DT * SUBSTEPS_DT;
        particle.x_prev = oldx;
        particle.y_prev = oldy;
        particle.ax = 0.0;
        particle.ay = 0.0;

        if (particle.x < XMin)
        {
            particle.x = XMin;
            // particle.x = XMin;

            // particle.vx = -DAMPING * particle.vx;
        }
        if (particle.x > XMax)
        {
            particle.x = XMax;
            // particle.vx = -DAMPING * particle.vx;
        }
        if (particle.y < YMin)
        {
            particle.y = YMin;
            // particle.vy = -DAMPING * particle.vy;
        }
        if (particle.y > YMax)
        {
            particle.y = YMax;
            // particle.vy = -DAMPING * particle.vy;
        }
    }
}
void init_particles(Particle *particles, sf::CircleShape *shapes, int count)
{
    const float center_x = (XMax + XMin) * 0.5f;
    const float center_y = (YMax + YMin) * 0.5f;
    const float span_x = (XMax - XMin) * 0.5f;
    const float span_y = (YMax - YMin) * 0.5f;
    const float span = (span_x > span_y ? span_x : span_y);
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        
        // particle.x = randf(XMin*MAGNIFICATION, XMax*MAGNIFICATION);
        // particle.y = randf(YMin*MAGNIFICATION, YMax*MAGNIFICATION);
        // Initialize velocity to be in stable circular motion. centripital == gravitational force
        // float angle = randf(0.0f, 2.0f * 3.14159265358979323846f);

        float x = center_x + randf(-OUTER_RADII * span_x, OUTER_RADII * span_x);
        float y = center_y + randf(-OUTER_RADII * span_y, OUTER_RADII * span_y);
        while ((NORM2(x - center_x, y - center_y) < SQ(INNER_RADII * span)) || NORM2(x - center_x, y - center_y) > SQ(OUTER_RADII * span))
        {
            x = center_x + randf(-OUTER_RADII * span_x, OUTER_RADII * span_x);
            y = center_y + randf(-OUTER_RADII * span_y, OUTER_RADII * span_y);
        }
        // const float radius = randf(INNER_RADII * (XMax - XMin) * 0.5f, OUTER_RADII * (XMax - XMin) * 0.5f); // Use half the width as radius
        particle.x = x;
        particle.y = y;
        particle.x_prev = particle.x; 
        particle.y_prev = particle.y; 
    }
    const float initial_iterations = 10;
    for(int i = 0; i < initial_iterations; i++)
    {
        float bounds[3];
        get_bounding_box(particles, count, bounds);
        tree_array_main = construct_trees(particles, tree_array_main, count, bounds);
        // print(tree_array_main[0].valid ? "Root node is valid." : "Root node is not valid.");
        const int ncount = current_available_index;
        computer_tree_coms(particles, tree_array_main, count, ncount);
        handle_collisions(particles, tree_array_main, count, ncount, COLLISION_HANDLING_PER_ITERATION);
        

    }
    for(int i =0; i<10; i++){
    step_particles(particles, count);
    }
    stop_particles(particles, count);
    
    for(int i = 0; i < count; i++){
        Particle &particle = particles[i];
        sf::CircleShape &shape = shapes[i];

        // particle.x = ((XMax + XMin) / 2.0f + radius * cosf(angle));
        // particle.y = ((YMax + YMin) / 2.0f + radius * sinf(angle));
        // particle.x = (XMax + XMin) / 2.0f + RAD * cosf(angle);
        // particle.y = (YMax + YMin) / 2.0f + RAD * sinf(angle);
        float angle = atan2f(particle.y - center_y, particle.x - center_x);                                                                           // Angle from the center to the particle
        float speed = CLOUD_VEL_MULTIPLIER * sqrtf(G * CENTER_MASS * CENTER_MASS_COUNT / sqrtf(NORM2(particle.x - center_x, particle.y - center_y))); // Speed for stable circular motion
        
        const double vx = speed * sin(angle);
        const double vy = -speed * cos(angle);
        particle.x_prev = particle.x - vx * SUBSTEPS_DT; // Previous position for Verlet integration
        particle.y_prev = particle.y - vy * SUBSTEPS_DT;

        // particle.x_prev = particle.x; // Previous position for Verlet integration
        // particle.y_prev = particle.y;

        particle.ax = 0.0;
        particle.ay = 0.0;

        particle.dx = 0.0;
        particle.dy = 0.0;
        // particle.dvx = 0.0;
        // particle.dvy = 0.0;
        particle.rad = RAD + randf(-dRad, dRad);

        // particle.mass = RAD_MASS_CONSTANT * particle.rad * particle.rad;
        particle.mass = MASS;
        shape.setRadius(1.0);
        shape.setOrigin(particle.rad, particle.rad);
        shape.setPosition(particle.x, particle.y);
#ifdef RAD_LUMINOSITY
        shape.setFillColor(sf::Color(255, 255, 255, (int)round(SQ(particles[i].rad) * RAD_LUMINOSITY)));
#else
        shape.setFillColor(sf::Color::White);
#endif
    }
    particles[0].x = D_CENTER + ((XMax + XMin) / 2.0f);
    particles[0].y = ((YMax + YMin) / 2.0f);
    particles[1].x = -D_CENTER + ((XMax + XMin) / 2.0f);
    particles[1].y = ((YMax + YMin) / 2.0f);
    float vy = sqrtf(G * CENTER_MASS / (4* D_CENTER));
    vy*= -1;
    particles[1].y_prev = particles[1].y - vy * SUBSTEPS_DT; // Previous position for Verlet integration
    particles[1].x_prev = particles[1].x; // Previous position for Verlet integration

    // particles[0].vy = -particles[1].vy;
    vy = -vy; // Reverse the velocity for the second particle
    particles[0].y_prev = particles[0].y - vy * SUBSTEPS_DT; // Previous position for Verlet integration
    particles[0].x_prev = particles[0].x; // Previous position for Verlet integration
    for (int i = 0; i < 2; i++)
    {
        particles[i].rad = CENTER_RAD;
        particles[i].mass = CENTER_MASS;
        // particles[i].vx = 0.0f;
        // particles[i].x_prev = particles[i].x;
        particles[i].ax = 0.0f;
        particles[i].ay = 0.0f;
        shapes[i].setRadius(particles[i].rad);
        shapes[i].setOrigin(particles[i].rad, particles[i].rad);
        shapes[i].setPosition(particles[i].x, particles[i].y);

        shapes[i].setFillColor(sf::Color::Red);
    }

    std::cout << "Particles initialized." << std::endl;
}

void init_velocities(Particle *particles, int count)
{
        const float center_x = (XMax + XMin) * 0.5f;
    const float center_y = (YMax + YMin) * 0.5f;
    const float span_x = (XMax - XMin) * 0.5f;
    const float span_y = (YMax - YMin) * 0.5f;
    // Initialize velocities for particles
    for(int i = 0; i < count; i++){
        Particle &particle = particles[i];

        // particle.x = ((XMax + XMin) / 2.0f + radius * cosf(angle));
        // particle.y = ((YMax + YMin) / 2.0f + radius * sinf(angle));
        // particle.x = (XMax + XMin) / 2.0f + RAD * cosf(angle);
        // particle.y = (YMax + YMin) / 2.0f + RAD * sinf(angle);
        float angle = atan2f(particle.y - center_y, particle.x - center_x);                                                                           // Angle from the center to the particle
        float speed = CLOUD_VEL_MULTIPLIER * sqrtf(G * CENTER_MASS * CENTER_MASS_COUNT / sqrtf(NORM2(particle.x - center_x, particle.y - center_y))); // Speed for stable circular motion
        
        const double vx = speed * sin(angle);
        const double vy = -speed * cos(angle);
        particle.x_prev = particle.x - vx * SUBSTEPS_DT; // Previous position for Verlet integration
        particle.y_prev = particle.y - vy * SUBSTEPS_DT;

        // particle.x_prev = particle.x; // Previous position for Verlet integration
        // particle.y_prev = particle.y;

        particle.ax = 0.0;
        particle.ay = 0.0;

        particle.dx = 0.0;
        particle.dy = 0.0;

    }
    particles[0].x = D_CENTER + ((XMax + XMin) / 2.0f);
    particles[0].y = ((YMax + YMin) / 2.0f);
    particles[1].x = -D_CENTER + ((XMax + XMin) / 2.0f);
    particles[1].y = ((YMax + YMin) / 2.0f);
    float vy = sqrtf(G * CENTER_MASS / (4* D_CENTER));
    vy*= -1;
    particles[1].y_prev = particles[1].y - vy * SUBSTEPS_DT; // Previous position for Verlet integration
    particles[1].x_prev = particles[1].x; // Previous position for Verlet integration

    // particles[0].vy = -particles[1].vy;
    vy = -vy; // Reverse the velocity for the second particle
    particles[0].y_prev = particles[0].y - vy * SUBSTEPS_DT; // Previous position for Verlet integration
    particles[0].x_prev = particles[0].x; // Previous position for Verlet integration
    for (int i = 0; i < 2; i++)
    {
        particles[i].rad = CENTER_RAD;
        particles[i].mass = CENTER_MASS;
        // particles[i].vx = 0.0f;
        // particles[i].x_prev = particles[i].x;
        particles[i].ax = 0.0f;
        particles[i].ay = 0.0f;
       
    }
}


// void step_gravity(Particle *particles, int count)
// {
//     // Verlet integration for gravity between particles
//     for (int i = 0; i < count; i++)
//     {
//         Particle &particle = particles[i];
//         float fx = 0.0f;
//         float fy = 0.0f;

//         for (int j = i; j < count; j++)
//         {
//             if (i == j)
//                 continue; // Skip self

//             Particle &other = particles[j];
//             float dx = (other.x - particle.x);
//             float dy = (other.y - particle.y);
//             float dist_sq = dx * dx + dy * dy;

//             if (dist_sq < 1e-2f)
//                 continue; // Avoid division by zero

//             const float dist = sqrtf(dist_sq);
//             float force = G * (particle.mass * other.mass) / dist_sq;
//             force = force > MAX_FORCE ? 0 : force; // Limit force to MAX_FORCE
//             particle.ax += force * (dx / dist) / particle.mass;
//             particle.ay += force * (dy / dist) / particle.mass;
//             other.ax -= force * (dx / dist) / other.mass;
//             other.ay -= force * (dy / dist) / other.mass;
//         }
//     }
//     // Update positions and velocities
//     for (int i = 0; i < count; i++)
//     {
//         Particle &particle = particles[i];
//         particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
//         particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
//         // if(i==100) std::cout << "Particle " << i << " position: (" << particle.x/MAGNIFICATION << ", " << particle.y/MAGNIFICATION << ")" << std::endl;
//         particle.vx += particle.ax * SUBSTEPS_DT;
//         particle.vy += particle.ay * SUBSTEPS_DT;
//         particle.x += particle.vx * SUBSTEPS_DT * 0.5f;
//         particle.y += particle.vy * SUBSTEPS_DT * 0.5f;
//         particle.ax = 0.0f;
//         particle.ay = 0.0f;
//     }
// }

void sync_shapes(Particle *particles, sf::CircleShape *shapes, int count)
{
    for (int i = 0; i < count; i++)
    {
        Particle &particle = particles[i];
        sf::CircleShape &shape = shapes[i];
        shape.setPosition(particle.x, particle.y);
    }
}


int count = 0;
void perform_barnes_hut(Particle *particles, Node *tree_array, const int pcount, const float theta)
{
    float bounds[3];
    print("Performing Barnes-Hut algorithm...");
    TIMEIT(get_bounding_box(particles, pcount, bounds), "Getting bounding box");
    print("Constructing trees...");

    TIMEIT(tree_array_main = construct_trees(particles, tree_array_main, pcount, bounds), "Constructing trees");
    print("Tree checking...");
    print(tree_array_main[0].valid ? "Root node is valid." : "Root node is not valid.");
    const int ncount = current_available_index;
    print("Trees constructed successfully. Number of nodes: " + std::to_string(ncount));
    print("Computing center of mass for nodes...");
    TIMEIT(computer_tree_coms(particles, tree_array_main, pcount, ncount), "Computing center of mass");
    print("Center of mass computed successfully.");
    print("Calculating accelerations using Barnes-Hut algorithm...");
    TIMEIT(calculate_accelerations_barnes_hut(particles, tree_array_main, pcount, ncount, theta), "Calculating accelerations");
    TIMEIT(handle_collisions(particles, tree_array_main, pcount, ncount, COLLISION_HANDLING_PER_ITERATION), "Handling collisions");
    print("Accelerations calculated successfully.");
    if(count==10){
        stop_particles(particles, pcount);
        init_velocities(particles, pcount);
    }
    count++;
}

int main()
{
    // Initialize random seed
    srand(42);
    int frame_count = 10000;

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
            // print("Performing Barnes-Hut algorithm...");
            perform_barnes_hut(particles_main, tree_array_main, PARTICLE_COUNT, THETA);
            TIMEIT(step_particles(particles_main, PARTICLE_COUNT), "Stepping particles");
        }
        sync_shapes(particles_main, pshape_main, PARTICLE_COUNT);

        for (int i = 0; i < PARTICLE_COUNT; i++)
        {
            window.draw(pshape_main[i]);
        }
        sf::Texture texture;
        texture.create(window.getSize().x, window.getSize().y);
        texture.update(window);
        sf::Image screenshot = texture.copyToImage();
        std::string filename = "frames/particles-" + std::to_string(PARTICLE_COUNT) + "_substeps-dt-" + std::to_string(SUBSTEPS_DT) + "_SUBSTEPS-" + std::to_string(SUBSTEPS) + "_theta-" + std::to_string(THETA) + "_frame-" + std::to_string(frame_count) + ".png";
        // screenshot.saveToFile(filename);
        frame_count++;
        window.display();
    }
    return 0;
}