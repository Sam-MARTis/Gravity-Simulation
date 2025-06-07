
#include <iostream>
#include <stdlib.h>


void print(auto message){
    std::cout << message << std::endl;
}
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
    int leaf1idx;
    int leaf2idx;
    int leaf3idx;
    int leaf4idx;
    float centerX;
    float centerY;
    float sideLength;
    
    Node() : valid(false), leaf1idx(-1), leaf2idx(-1), leaf3idx(-1), leaf4idx(-1), centerX(0.0f), centerY(0.0f), sideLength(0.0f) {}
    Node(bool v, int l1, int l2, int l3, int l4, float cx, float cy, float sl) 
        : valid(v), leaf1idx(l1), leaf2idx(l2), leaf3idx(l3), leaf4idx(l4), centerX(cx), centerY(cy), sideLength(sl) {}
    Node(const Node& other)
        : valid(other.valid), leaf1idx(other.leaf1idx), leaf2idx(other.leaf2idx), leaf3idx(other.leaf3idx), leaf4idx(other.leaf4idx),
          centerX(other.centerX), centerY(other.centerY), sideLength(other.sideLength) {}
};


int main(){

    Node n;
    print(n.valid);
    print(sizeof(n));
    print(sizeof(Node));
    print(sizeof(bool) + sizeof(int) * 4 + sizeof(float) * 3);

    print(*(int*)(&n.valid + 4*sizeof(bool)));
    print(*&n.leaf1idx);
    print(&n.leaf2idx);
    print(&n.leaf3idx);

    
    return 0;
}