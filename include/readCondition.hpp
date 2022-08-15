#ifndef READ_CONDITION_HPP
#define READ_CONDITION_HPP

#include <fstream>
#include <unordered_map>
#include <vector>

#include <spdlog/spdlog.h>

#include "info.hpp"

namespace lhfea{

bool readCondition(std::string path, realn &E, realn &nu,
            std::unordered_map<int,std::vector<realn>> &bound)
{
    std::ifstream in(path);
    in >> E >> nu;
    try{
        while(!in.fail()){
            int idx; in >> idx;
            realn tmp;
            std::vector<realn> list;
            for(int i=0; i<9; i++){
                in >> tmp;
                list.push_back(tmp);
            }
            bound.emplace(idx, list);
        }
    }
    catch(...){
        spdlog::error("The input file does not meet the format requirements.");
        exit(-1);
    }
    in.close();
    return true;
}

} // namespace lhfea


#endif