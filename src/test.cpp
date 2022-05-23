#include <iostream>
#include <string>
int main(){
    int i = 5;
    std::string str_ = "wall_temp_" + std::to_string(i);
    std::cout<<str_<<"\n";
}