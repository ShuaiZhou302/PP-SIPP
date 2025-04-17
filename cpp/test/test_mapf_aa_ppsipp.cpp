/*******************************************
 * Author: Shuai Zhou.
 * All Rights Reserved.
 *******************************************/

#include "mapfaa_ppsipp.hpp"
#include "debug.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip> // Include this header for std::fixed and std::setprecision
#include <unordered_map>

int TestPP_SIPP();
int edge_hash(long a, long b);
int Toyexample();

//int main(){
//    TestPP_SIPP();
//    return 0;
//};


int TestPP_SIPP(){
    std::cout << "####### test_pibt.cpp - TestLsrp() Begin #######" << std::endl;
    // ------------static environment, 3x3, obstacle at (y=0,x=0),(y=0.x=1),(y=0,x=2),(y=1,x=0) Three robots---------
    //   x:0  1  2
    // y ---------------
    // 0 | *  G2 *  |
    // 1 | R1 *  *  |
    // 2 | *  G1 *  |
    //   ---------------
    //  Note: *=free, #=obstacle, R1,R2,R3 = robot start, G1,G2,G3 = robot goals
    // -------------------------------------------------------------------------------
    raplab::Grid2d g;
    int r=3,c=3;
    std::vector<std::vector<double>> occupancy_grid;
    occupancy_grid.resize(r);
    for (int i = 0; i < r; i++){
        occupancy_grid[i].resize(c, 0);
    }
//    occupancy_grid[1][0] = 1;
//    occupancy_grid[1][2] = 1;
//    occupancy_grid[2][0] = 1;
//    occupancy_grid[2][2] = 1;
//    occupancy_grid[3][0] = 1;
//    occupancy_grid[3][2] = 1;
    //occupancy_grid[1][0] = 1;
    g.SetOccuGridPtr(&occupancy_grid);
    // not using it now
    double time_limit = 3600;
//    std::vector<long> starts({3,7,5}); // node id = y*NumX + x; e.g. (y=1)*4 + (x=0) = 4
    std::vector<long> starts({7,3});
    // the length of starts shows how many agents are there.
//    std::vector<long> goals({7,5,3});
    std::vector<long> goals({1,7});
    raplab::PP_SIPP planner;
//    std::vector<double> dura({1,2,3});
    std::vector<double> duration({2,1});
    planner.SetGraphPtr(&g);
    //std::vector<double> duration(starts.size(),1);
    //std::vector<double> duration = {1,2};
    planner.Solve(starts, goals, time_limit, 1.0);
    planner.set_duration(duration);
    auto result = planner.run_PP();
    if (result) {
        std::cout << "Path found!" << std::endl;
        std::cout << "Cost "<<planner.get_sum_of_cost()<<std::endl;
    } else {
        std::cout << "No path found." << std::endl;
    }
    std::cout << "####### test_LSRP.cpp - TestLSRP() End #######" << std::endl;
    return 0;
}

int Toyexample() {
    std::cout << "####### Toy-example Begin #######" << std::endl;
    /*
    3x3 Grid graph：
    (0,0) (0,1) (0,2)                           0  1  2
    (1,0) (1,1) (1,2)  and corresponding ids:   3  4  5
    (2,0) (2,1) (2,2)                           6  7  8
    illustration：2D coordinate (x, y) can be transformed to id：ID = y * width + x
    eg： (0,0) -> 0, (2,2) -> 8

    Obstacles position (X) ：
      0  1  2
      X  4  5
      6  7  8

    start and goals (S) (G) ：
      S1  1  G2
      X   4  5
      S2  7  G1

    time for each agent to move from one cell to another:
    agent 1 duration = 1, agent 2 duration = 0.5

    specific edge duration for specific agent (optional)：
    agent 1： duration = 0.5 at edge (0,1) to (0,2)
    */

    raplab::Grid2d g;
    std::vector<std::vector<double>> occupancy_grid(3, std::vector<double>(3, 0));
    // set
    occupancy_grid[1][0] = 1; // obstacle 1
    g.SetOccuGridPtr(&occupancy_grid);
    // set starts and goals
    std::vector<long> starts = {0, 6}; // (0,0) and (2,0)
    std::vector<long> goals = {8, 2};  // (2,2) and (0,2)
    // set duration for each agent
    std::vector<double> duration = {1, 0.5}; //
    // set specific edge duration for specific agent if required (optional)
    std::unordered_map<int, std::unordered_map<int, double>> edge_duration = {
            {1, {{edge_hash(1, 2), 0.5}}}
    }; // agent 1 duration = 0.5 at edge (0,1) to (0,2)
    // set time limit
    double time_limit = 30;
    // initialize
    raplab::PP_SIPP planner;
    planner.SetGraphPtr(&g);
    // search
    planner.Solve(starts, goals, time_limit, 1.0);  // eps is not used
    planner.set_duration(duration);
    planner.set_edge_cost(edge_duration); // optional
    auto succ = planner.run_PP();
    // get the result
    auto soc = planner.get_sum_of_cost();
    auto makespan = planner.get_makespan();
    auto runtime = planner.get_runtime();
    // output the result
    if (succ) {
        std::cout << "Solution found: true" << std::endl;
        std::cout << "Runtime: " << std::fixed << std::setprecision(3) << runtime << std::endl;
        std::cout << "Makespan: " << std::fixed << std::setprecision(2) << makespan << std::endl;
        std::cout << "Soc: " << std::fixed << std::setprecision(2) << soc << std::endl;
    } else {
        std::cout << "Solution found: false" << std::endl;
    }
    std::cout << "####### Toy-example End #######" << std::endl;
    return 1;
}


int edge_hash(long a, long b)
{
    //use cantor hash to match the method inside lsrp
    if (a > b)
    {
        return(a + b) * (a + b + 1) / 2 + b;
    } else {
        return(a + b) * (b + a + 1) / 2 + a;
    }
}