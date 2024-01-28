#include <iostream>
using namespace std;

int main() {
    int count {};
    string name {};
    cout << "Enter name: ";
    cin >> name;
    cout << "Enter positive number: ";
    cin >> count;
    
    string tag {};

    if (count == 1) {
        tag = "st";
    }
    else if (count == 2) {
        tag = "nd";
    }
    else if (count == 3) {
        tag = "rd";
    }
    else {
        tag = "th";
    }
    cout << "Hello " << name << "! This is your " << count << tag << " program in C++ from Ubuntu" << endl;

    return 0;
}