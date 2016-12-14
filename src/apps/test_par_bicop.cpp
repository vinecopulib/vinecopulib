//
// Created by Vatter Thibault on 14/12/16.
//


#include <common/include/par_bicop.h>
#include <iostream>
#include <vector>

using namespace std;

int main(int __unused argc, char __unused *argv[]) {


    vector<int> allfams = {0,1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,23,24,26,27,28,29,30,23,33,34,36,37,38,39,
                                 40,104,114,124,134,204,214,224,234};
    vector<int> tawns = {33,34,35,36,37,38,39,40};
    vector<int> onepar = {1,3,4,5,6,11,12,13,18,19,20,25,26,27};
    vector<int> twopar = {2,7,8,9,10,14,15,16,17,21,22,23,24,28,29,30,31,32,33,34,35,36,37,38,39};
    vector<int> negfams = {1,2,5,23,24,26,27,28,29,30,33,34,36,37,38,39,40,124,134,224,234};
    vector<int> posfams = {1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,104,114,204,214};

    // Iterate and print values of vector
    cout << "All available families: \n";
    for(int n : allfams) {
        cout << n << " ";
    }
    cout << "\n One-parameter families: \n";
    for(int n : onepar) {
        cout << allfams[n] << " ";
    }
    cout << "\n Two-parameters families: \n";
    for(int n : twopar) {
        cout << allfams[n] << " ";
    }
    cout << "\n Families for positive dependence: \n";
    for(int n : posfams) {
        cout << n << " ";
    }
    cout << "\n Families for negative dependence: \n";
    for(int n : negfams) {
        cout << n << " ";
    }
    cout << "\n";
    return 0;
}
