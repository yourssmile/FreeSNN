// #ifndef  _POISSON_PROCESS_H
// #define  _POISSON_PROCESS_H

// #include <array>
// #include <fstream>
// #include <iostream>
// #include <random>
// #include <vector>

// class PoissonProcess {
// public:
//     void run() {
//         std::random_device rd;
//         std::mt19937 gen(rd());
//         std::uniform_real_distribution<double> dis(0.0, 1.0);
    
//         double cur = 0.0;
//         while(cur < total) {
//             double r = dis(gen);
//             double dt = -(1.0 / lamda) * log(r);
//             cur = cur + dt;
//             count++;
//             if(cur < total) {
//                 occur.push_back(cur);
//             }
//         }
        
//         std::cout << "count = " << count << std::endl;
//         std::cout << "total = " << total << std::endl;
//     }

//     void print() {
//         for(auto& e : occur) {
//             std::cout << e << std::endl;
//         }
//     }

//     void verify() {
//         std::array<std::vector<double>, num> slices;

//         for(auto&e : occur) {
//             int index = floor(e * num / total);
//             slices[index].push_back(e);
//         }
        
//         double expect = lamda * total / num;
//         int range = ceil(3*expect);  // experience
//         std::vector<int> data(range, 0);
        
//         int few = 0;
//         for(int i = 0; i < num; ++i) {
//             if(slices[i].size() >= range) { // out of bound
//                 few++;
//             }
//             else {
//                 data[slices[i].size()]++;
//             }
//         }
//         std::cout << "few = " << few << std::endl;

//         std::ofstream of;
//         of.open("result.txt");
//         double factor = 1.0;
//         for (int i = 0; i < range; i++) {
//             if(i != 0) {
//                 factor = factor * i;
//             }
    
//             auto theory_out = exp(-1 * expect) * pow(expect, i) / factor;
//             auto sim_out = (double)data[i] / num;
//             of << i << " " << sim_out << " " << theory_out << std::endl;
//         }
//         of.close();
//     }

// private:
//     static const int num = 2500;  // slice num
//     const double total = 1000.0;  // total time
//     const double lamda = 20.0;   // λ

//     int count;
//     std::vector<double> occur;
// };

// #endif // _POISSON_PROCESS_H