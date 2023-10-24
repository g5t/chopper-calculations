#include <iostream>
#include "choppers.h"
int main() {
    for (const auto & [name, value] : mcstas_calculation(7.0, 0., 0.0002)){
        std::cout << name << " " << value << std::endl;
    }
    return 0;
}
