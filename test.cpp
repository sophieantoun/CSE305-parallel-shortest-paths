#include <omp.h>
#include <iostream>

int main() {
    #pragma omp parallel
    {
        std::cout << "Hello from thread " << omp_get_thread_num() << std::endl;
    }
    return 0;
}
