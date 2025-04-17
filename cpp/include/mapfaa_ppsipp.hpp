/*******************************************
 * Author: Shuai Zhou.
 * All Rights Reserved.
 *******************************************/
#ifndef LSRP_LNS_MAPFAA_PPSIPP_HPP
#define LSRP_LNS_MAPFAA_PPSIPP_HPP

#include "mapfaa_util.hpp"
#include <utility>
#include <vector>
#include <tuple>
#include <queue>
#include <unordered_map>
#include <random>
#include <optional>
#include <algorithm>
#include <limits>
#include <unordered_set>
#include <iostream>
#include <string>
#include <cassert>
#include <memory>
#include <algorithm>
#include <set>
#include <map>

namespace raplab{

    struct Sipp_Node{
        long _v;
        double _t;
        std::pair<double, double> _interval;
        double _f;
        Sipp_Node* _parent = nullptr;

        Sipp_Node() = default;
        Sipp_Node(long v,double t, std::pair<double, double> interval,double f, Sipp_Node* parent = nullptr)
                : _v(v), _t(t), _interval(std::move(interval)), _f(f), _parent(parent) {}

        bool operator==(const Sipp_Node& other) const {
            return _v == other._v &&  _interval.first == other._interval.first
                   && _interval.second == other._interval.second;
        }

        struct Hash {
            std::size_t operator()(const Sipp_Node& node) const {
                std::size_t h1 = std::hash<long>()(node._v);               // Hash for vertex
                std::size_t h2 = std::hash<double>()(node._interval.first); // Hash for interval start
                std::size_t h3 = std::hash<double>()(node._interval.second); // Hash for interval end
                return h1 ^ (h2 << 1) ^ (h3 << 2); // Combine the hashes
            }
        };
    };

    struct CompareF {
        bool operator()(const Sipp_Node* a, const Sipp_Node* b) const {
            if (a->_f == b->_f) {
                return a->_t < b->_t; // Sort by _t if _f is equal
            }
            return a->_f < b->_f; // Sort by _f primarily
        }
    };

    struct PairHash {
        std::size_t operator()(const std::pair<long, std::pair<double, double>>& p) const {
            std::size_t h1 = std::hash<long>()(p.first);
            std::size_t h2 = std::hash<double>()(p.second.first);
            std::size_t h3 = std::hash<double>()(p.second.second);
            return h1 ^ (h2 << 1) ^ (h3 << 2); // Combine the hashes
        }
    };
    class PP_SIPP : public MAPFAAPlanner {
    public:
        PP_SIPP();
        virtual ~PP_SIPP();

        virtual CostVec GetPlanCost(long nid=-1) override ;
        virtual TimePathSet GetPlan(long nid=-1) override ;
        virtual std::unordered_map<std::string, double> GetStats() override ;
        virtual int Solve(std::vector<long>& starts, std::vector<long>& goals, double time_limit, double eps) override;
        void set_duration(std::vector<double> duration);
        void set_edge_cost(std::unordered_map<int,std::unordered_map<int,double>> edge_cost) {this->edge_cost = edge_cost;}
        double get_h(int agent, long vertex);
        double get_D(int agent, long v_departure,long v_arrive);
        bool run_PP();
        bool run_SIPP(int agent);

        std::vector<std::unordered_map<long,double>> generate_distable();

        std::unordered_map<long,double> generate_single_dis_table(int agent);

        std::vector<std::pair<double, double>> get_safe_intervals(long vertex);

        std::vector<Sipp_Node> get_successors(Sipp_Node* n, int agent);

        double get_sum_of_cost() {return _cost;}
        double get_makespan() { return  _makespan;}
        double get_runtime() { return _runtime; }
        std::vector<std::vector<std::tuple<long, long, double, double>>> get_all_paths() { return _all_paths; }
        bool Timeout() {
            auto curr = std::chrono::steady_clock::now();
            return std::chrono::duration_cast<std::chrono::duration<double>>(curr - _start_time).count() >= _time_limit;
        }

        void remove_node_from_open(Sipp_Node* existing_node) {
            for (auto it = Open.begin(); it != Open.end(); ++it) {
                if (**it == *existing_node) { // compare
//                    delete *it;
                    it = Open.erase(it); // remove the node from Open
                }
            }
        }

        int edge_hash(long a, long b) const{
            if (a > b)
            {
                return(a + b) * (a + b + 1) / 2 + b;
            }
            return(a + b) * (b + a + 1) / 2 + a;
        }


        void reconstruct_Path(int agent, Sipp_Node* n);

        void add_to_unsafe_intervals(long vertex, std::pair<double, double> interval) {
            if (_unsafe_intervals.find(vertex) == _unsafe_intervals.end()) {
                _unsafe_intervals[vertex] = std::set<std::pair<double, double>, IntervalCompare>();
            }
            _unsafe_intervals[vertex].insert(interval);
        }

        void initialize_unsafe_intervals(const std::vector<std::vector<std::tuple<long, long, double, double>>>& input_path) {
            _unsafe_intervals.clear();
            for (const auto& path : input_path) {
                for (const auto& move : path) {
                    auto p = std::get<0>(move);
                    auto v = std::get<1>(move);
                    auto tp = std::get<2>(move);
                    auto tv = std::get<3>(move);
                    add_to_unsafe_intervals(p, {tp, tv});
                    if (move == *std::prev(path.end())) {
                        // last move
                        add_to_unsafe_intervals(v, {tp, std::numeric_limits<double>::infinity()});
                    } else {
                        add_to_unsafe_intervals(v, {tp, tv});
                    }
                }
            }
        }




    private:
        TimePathSet _paths;
        std::vector<std::unordered_map<long,double>> _dis_table;
        std::vector<std::vector<std::tuple<long, long, double, double>>> _all_paths;
        std::vector<double> _duration;
        std::unordered_map<int,std::unordered_map<int,double>> edge_cost;
        double _time_limit;
        double _runtime;
        std::chrono::steady_clock::time_point _start_time;
        double _cost;
        double _makespan;

        //std::priority_queue<Sipp_Node*, std::vector<Sipp_Node*>, CompareF> Open;
        std::multiset<Sipp_Node*, CompareF> Open;
        std::unordered_set<Sipp_Node, Sipp_Node::Hash> Closed;
        std::unordered_map<std::pair<long, std::pair<double, double>>, Sipp_Node, PairHash> Visited;

        struct IntervalCompare {
            bool operator()(const std::pair<double, double>& a, const std::pair<double, double>& b) const {
                return a.first < b.first;
            }
        };
        std::unordered_map<long, std::set<std::pair<double, double>, IntervalCompare>> _unsafe_intervals;
        //std::unordered_map<long, std::vector<std::pair<double, double>>> _safe_intervals;
    };

}
#endif //LSRP_LNS_MAPFAA_PPSIPP_HPP
