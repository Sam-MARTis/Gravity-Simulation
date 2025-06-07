
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
   
};


int main(){

    // Node n;
    // print(n.valid);
    // print(sizeof(n));
    // print(sizeof(Node));
    // print(sizeof(bool) + sizeof(int) * 4 + sizeof(float) * 3);

    // print(*(int*)(&n.valid + 4*sizeof(bool)));
    // print(*&n.leaf1idx);
    // print(&n.leaf2idx);
    // print(&n.leaf3idx);

    Node *tree = new Node[100];
    for(int i = 0; i < 100; i++){
        print(tree[i].valid);
    }


    
    return 0;
}